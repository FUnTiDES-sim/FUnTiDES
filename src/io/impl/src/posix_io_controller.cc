#include "posix_io_controller.h"

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace funtides::io {
namespace {

constexpr char kMagicSnap[8] = {'F', 'U', 'N', 'T', 'S', 'N', 'A', 'P'};
constexpr std::uint32_t kFormatVersion = 1;
constexpr std::size_t kMaxDims = 4;

/// Fixed-size, naturally aligned. Written verbatim at the start of each file.
struct FileHeader {
  char magic[8];
  std::uint32_t version;
  std::uint32_t ndim;
  std::uint64_t dims[kMaxDims];         ///< Local extents of the payload.
  std::uint64_t global_dims[kMaxDims];  ///< Global shape, 0 if not distributed.
  std::uint64_t offsets[kMaxDims];      ///< This rank's offset in the global shape.
  std::uint64_t snapshot_index;         ///< Ordinal of this snapshot in the run.
  std::uint32_t scalar_bytes;
  std::uint32_t reserved;
  std::uint64_t nelem;  ///< Number of scalars in the payload.
};
static_assert(sizeof(FileHeader) == 136, "unexpected padding in FileHeader");

/// `shot_id` becomes a path component, so anything that could escape the output
/// directory is rejected rather than sanitized silently.
void validateShotId(const std::string& shot_id) {
  for (const char c : shot_id) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
    if (!ok) {
      throw std::invalid_argument("funtides::io: shot_id must contain only alphanumerics, '_' or '-', got \"" +
                                  shot_id + "\"");
    }
  }
}

class File {
 public:
  File(const std::string& path, const char* mode) : path_(path) {
    fp_ = std::fopen(path.c_str(), mode);
    if (fp_ == nullptr) {
      throw std::runtime_error("funtides::io: cannot open " + path + ": " + std::strerror(errno));
    }
  }
  ~File() {
    if (fp_ != nullptr) std::fclose(fp_);
  }

  File(const File&) = delete;
  File& operator=(const File&) = delete;

  void write(const void* data, std::size_t bytes) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp_) != bytes) {
      throw std::runtime_error("funtides::io: short write on " + path_);
    }
  }

  void read(void* data, std::size_t bytes) {
    if (bytes == 0) return;
    if (std::fread(data, 1, bytes, fp_) != bytes) {
      throw std::runtime_error("funtides::io: short read on " + path_ + " (truncated file?)");
    }
  }

  void closeChecked() {
    if (fp_ == nullptr) return;
    const int rc = std::fclose(fp_);
    fp_ = nullptr;
    if (rc != 0) {
      throw std::runtime_error("funtides::io: error closing " + path_);
    }
  }

 private:
  std::FILE* fp_{nullptr};
  std::string path_;
};

FileHeader makeHeader(const char (&magic)[8]) {
  FileHeader h{};
  std::memcpy(h.magic, magic, sizeof(h.magic));
  h.version = kFormatVersion;
  h.scalar_bytes = static_cast<std::uint32_t>(sizeof(float));
  return h;
}

void checkHeader(const FileHeader& h, const char (&magic)[8], const std::string& path) {
  if (std::memcmp(h.magic, magic, sizeof(h.magic)) != 0) {
    throw std::runtime_error("funtides::io: " + path + " is not a FUnTiDES file");
  }
  if (h.version != kFormatVersion) {
    throw std::runtime_error("funtides::io: " + path + " has format version " + std::to_string(h.version) +
                             ", expected " + std::to_string(kFormatVersion));
  }
  if (h.scalar_bytes != sizeof(float)) {
    throw std::runtime_error("funtides::io: " + path + " was written with a different scalar size");
  }
}

}  // namespace

PosixIOController::PosixIOController(OpenMode mode, const IOConfig& config) : IOControllerBase(config), mode_(mode) {
  validateShotId(config.shot_id);

  if (mode_ == OpenMode::kWrite) {
    // create_directories() creates parents too, so the shot subdirectory and
    // the output directory itself both land here.
    // NOTE: every rank races on this call. error_code keeps an already-existing
    // directory from throwing, which is enough on a local filesystem; a
    // parallel filesystem would want rank 0 to create it and a barrier after,
    // which needs the communicator DistributedContext does not carry yet.
    const std::string dir = outputDir();
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
      throw std::runtime_error("funtides::io: cannot create " + dir + ": " + ec.message());
    }
  }
}

PosixIOController::~PosixIOController() {
  try {
    close();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "[funtides::io] error during close: %s\n", e.what());
  } catch (...) {
    std::fprintf(stderr, "[funtides::io] unknown error during close\n");
  }
}

std::string PosixIOController::outputDir() const {
  const std::string& shot = config().shot_id;
  return shot.empty() ? config().output_dir : config().output_dir + "/" + shot;
}

std::string PosixIOController::snapshotPath(std::size_t index) const {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "/%s_snap_%06zu_r%04d.bin", config().prefix.c_str(), index, config().ctx.rank);
  return outputDir() + buf;
}

void PosixIOController::requireMode(OpenMode expected, const char* what) const {
  if (mode_ != expected) {
    throw std::runtime_error(std::string("funtides::io: ") + what +
                             " called on a controller opened in the "
                             "opposite mode");
  }
}

void PosixIOController::writeSnapshot(const HostVectorReal& field) {
  requireMode(OpenMode::kWrite, "writeSnapshot");
  if (closed_) throw std::runtime_error("funtides::io: writeSnapshot after close");

  const IOConfig& cfg = config();

  FileHeader h = makeHeader(kMagicSnap);
  h.ndim = static_cast<std::uint32_t>(cfg.local_dims.size());
  if (h.ndim > kMaxDims) {
    throw std::runtime_error("funtides::io: local_dims has too many dimensions");
  }
  for (std::size_t d = 0; d < h.ndim; ++d) {
    h.dims[d] = cfg.local_dims[d];
    h.global_dims[d] = d < cfg.global_dims.size() ? cfg.global_dims[d] : 0;
    h.offsets[d] = d < cfg.start_offsets.size() ? cfg.start_offsets[d] : 0;
  }
  h.snapshot_index = next_index_;
  h.nelem = field.size();

  File f(snapshotPath(next_index_), "wb");
  f.write(&h, sizeof(h));
  f.write(field.data(), field.size() * sizeof(float));
  f.closeChecked();

  ++next_index_;
}

void PosixIOController::readSnapshot(const HostVectorReal& field, std::size_t index) {
  requireMode(OpenMode::kRead, "readSnapshot");

  const std::string path = snapshotPath(index);
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("funtides::io: no snapshot with index " + std::to_string(index) + " (" + path + ")");
  }

  File f(path, "rb");
  FileHeader h{};
  f.read(&h, sizeof(h));
  checkHeader(h, kMagicSnap, path);

  if (h.nelem != field.size()) {
    throw std::runtime_error("funtides::io: " + path + " holds " + std::to_string(h.nelem) + " values, the view has " +
                             std::to_string(field.size()));
  }
  f.read(field.data(), field.size() * sizeof(float));
}

void PosixIOController::flush() {
  // Each write opens, writes and closes its own file, so nothing is ever
  // pending at the library level. Data already left the process; whether it
  // reached the platters is up to the OS page cache, as everywhere else.
}

void PosixIOController::close() {
  if (closed_) return;
  closed_ = true;
  flush();
}

}  // namespace funtides::io
