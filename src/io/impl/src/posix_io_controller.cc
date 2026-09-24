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

/**
 * @brief On-disk header written verbatim (host byte order) at the start of each snapshot file.
 *
 * Fixed size, naturally aligned, followed by `nelem` scalars.
 */
struct FileHeader {
  char magic[8];
  std::uint32_t version;
  std::uint32_t ndim;
  std::uint64_t dims[kMaxDims];         ///< Local extents of the payload.
  std::uint64_t global_dims[kMaxDims];  ///< Global shape, 0 if not distributed.
  std::uint64_t offsets[kMaxDims];      ///< This rank's offset in the global shape.
  std::uint64_t snapshot_index;         ///< Ordinal of this snapshot in the run.
  std::uint32_t scalar_bytes;           ///< Size in bytes of one payload scalar.
  std::uint32_t reserved;
  std::uint64_t nelem;  ///< Number of scalars in the payload.
};
static_assert(sizeof(FileHeader) == 136, "unexpected padding in FileHeader");

/**
 * @brief Rejects a shot identifier that is not safe to use as a path component.
 *
 * `shot_id` becomes a directory name, so anything that could escape the output
 * directory is rejected rather than silently sanitized.
 *
 * @param[in] shot_id Candidate identifier; only alphanumerics, '_' and '-' are accepted.
 * @throws std::invalid_argument If a forbidden character is present.
 */
void validateShotId(const std::string& shot_id) {
  for (const char c : shot_id) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
    if (!ok) {
      throw std::invalid_argument("funtides::io: shot_id must contain only alphanumerics, '_' or '-', got \"" +
                                  shot_id + "\"");
    }
  }
}

/**
 * @brief RAII wrapper around a C `FILE*` that reports every failure by exception.
 *
 * Non-copyable. The destructor closes silently; call closeChecked() to detect
 * errors on close.
 */
class File {
 public:
  /// @throws std::runtime_error If the file cannot be opened.
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

  /// @throws std::runtime_error On a short write.
  void write(const void* data, std::size_t bytes) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp_) != bytes) {
      throw std::runtime_error("funtides::io: short write on " + path_);
    }
  }

  /// @throws std::runtime_error On a short read.
  void read(void* data, std::size_t bytes) {
    if (bytes == 0) return;
    if (std::fread(data, 1, bytes, fp_) != bytes) {
      throw std::runtime_error("funtides::io: short read on " + path_ + " (truncated file?)");
    }
  }

  /// @brief Closes the file and reports failure; no-op if already closed.
  /// @throws std::runtime_error If `fclose` fails.
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

/// @brief Returns a zero-initialized header with magic, version and scalar size filled in.
FileHeader makeHeader(const char (&magic)[8]) {
  FileHeader h{};
  std::memcpy(h.magic, magic, sizeof(h.magic));
  h.version = kFormatVersion;
  h.scalar_bytes = static_cast<std::uint32_t>(sizeof(float));
  return h;
}

/**
 * @brief Validates magic, format version and scalar size of a header read from `path`.
 * @throws std::runtime_error On any mismatch.
 */
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
    // create_directories() also creates the parents, so the output directory
    // and the shot subdirectory are both handled here.
    // TODO: every rank calls this concurrently; error_code keeps an already
    // existing directory from throwing, which is enough on a local filesystem.
    // A parallel filesystem needs rank 0 to create it followed by a barrier,
    // which requires a communicator that DistributedContext does not carry yet.
    const std::string dir = outputDir();
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
      throw std::runtime_error("funtides::io: cannot create " + dir + ": " + ec.message());
    }
  }
}

PosixIOController::~PosixIOController() {
  // A destructor must not throw: close errors are reported on stderr only.
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
  // Each write opens, writes and closes its own file, so nothing is pending
  // in this class. Durability beyond the OS page cache is not guaranteed.
}

void PosixIOController::close() {
  if (closed_) return;
  closed_ = true;
  flush();
}

}  // namespace funtides::io
