#ifndef FUNTIDES_IO_INCLUDE_IO_CONTROLLER_BASE_H_
#define FUNTIDES_IO_INCLUDE_IO_CONTROLLER_BASE_H_

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "data_type.h"
#include "distributed_ctx.h"

namespace funtides::io {

/** @brief Host mirror of a vectorReal: the buffer type of snapshot I/O. */
using HostVectorReal = vectorReal::host_mirror_type;
/**
 * @brief Host mirror of an arrayReal.
 *
 * No function of this module takes it; see docs/design-red-flags.md.
 */
using HostArrayReal = arrayReal::host_mirror_type;

/**
 * @brief Settings of one I/O controller: where files go and the shape of the
 * snapshots of this rank.
 *
 * Read and write controllers must be given the same `ctx.rank`, `output_dir`,
 * `prefix` and `shot_id` to address the same files.
 */
struct IOConfig {
  utils::DistributedContext ctx{};  ///< Rank of this process, part of file names.
  std::string output_dir{"."};      ///< Root output directory.
  std::string prefix{"funtides"};   ///< Prefix of every file name.

  /// Optional shot identifier. When set, output goes to <output_dir>/<shot_id>
  /// so that a multi-shot run keeps its files apart and no single directory
  /// accumulates millions of entries. Must contain only alphanumerics, '_' or
  /// '-': it becomes a path component.
  std::string shot_id{""};

  /// Extents of the global snapshot grid, one entry per axis.
  /// @todo VERIFY: axis order (x, y, z?) and unit (nodes per axis?).
  std::vector<std::size_t> global_dims;
  std::vector<std::size_t> start_offsets;  ///< Offset of this rank's chunk in global_dims, per axis.
  /// Extents of this rank's chunk, per axis; at most 4 axes.
  /// @todo VERIFY: must the size of every snapshot equal the product of
  /// local_dims? No backend checks it.
  std::vector<std::size_t> local_dims;

  /// @todo VERIFY: number of time steps or of snapshots? Read by no backend.
  std::size_t nt{0};
  /// @todo VERIFY: total number of receivers of the run? Read by no backend.
  std::size_t nb_receiver{0};

  /// Allows writeSnapshot() to return before the data is written. A backend
  /// may ignore it; the buffer contract of writeSnapshot() holds either way.
  bool async_snapshots{false};
};

/**
 * @brief Backend-agnostic interface to write and read back the wavefield
 * snapshots of one rank.
 *
 * Snapshots are numbered 0, 1, 2... in the order of the writeSnapshot() calls;
 * readSnapshot() takes that number. Create instances with makeIOController().
 *
 * HOST ONLY. Never capture an instance inside a Kokkos kernel.
 *
 * THREADING: not thread-safe. A single thread drives the controller.
 *
 * A snapshot is stored as the flat sequence of values of the view, in view
 * order; `local_dims` and `global_dims` are only recorded as metadata.
 * No solver or driver uses this interface yet; see docs/design-red-flags.md.
 */
class IOControllerBase {
 public:
  virtual ~IOControllerBase() = default;

  /**
   * @brief Writes the next snapshot of a field.
   *
   * Requires a controller opened with OpenMode::kWrite.
   * BUFFER CONTRACT: `field` must stay valid AND unmodified until the next
   * call to writeSnapshot(), flush(), or close(), because the backend may read
   * it after this call returns. Use two host mirrors alternately if the solver
   * needs to keep producing.
   * @todo VERIFY: which field (pressure, displacement components?) and in which
   * node ordering is a snapshot expected to hold?
   * @throws std::runtime_error on a controller opened for reading, after
   *         close(), or on an I/O error.
   */
  virtual void writeSnapshot(const HostVectorReal& field) = 0;

  /**
   * @brief Reads snapshot number `index` into a caller-allocated view.
   *
   * Requires a controller opened with OpenMode::kRead.
   * @param[out] field Destination; its size must equal the stored size.
   * @param[in] index Snapshot number, in writeSnapshot() call order.
   * @throws std::runtime_error if `index` is out of range, the extents of
   *         `field` do not match what is stored, or the controller was opened
   *         for writing.
   */
  virtual void readSnapshot(const HostVectorReal& field, std::size_t index) = 0;

  /**
   * @brief Blocks until every pending write has landed. Buffers handed to
   * writeSnapshot() are free afterwards.
   */
  virtual void flush() = 0;

  /**
   * @brief Flushes, then closes the files. Idempotent.
   *
   * Propagates errors, so call it explicitly: the destructor can only swallow
   * them. COLLECTIVE: every rank must call it, in the same order.
   */
  virtual void close() = 0;

 protected:
  /** @brief Stores a copy of `config`. */
  explicit IOControllerBase(const IOConfig& config) : config_(config) {}

  /** @brief Returns the settings given at construction. */
  const IOConfig& config() const { return config_; }

 private:
  IOConfig config_;
};

/**
 * @brief Storage back-end selected by makeIOController().
 *
 * kAdios2 has no implementation: makeIOController() throws for it; see
 * docs/design-red-flags.md.
 */
enum class BackendKind { kPosix, kAdios2 };
/** @brief Direction of a controller: it either writes or reads snapshots. */
enum class OpenMode { kWrite, kRead };

/**
 * @brief Creates an I/O controller.
 *
 * In OpenMode::kWrite the output directory, including the `shot_id`
 * subdirectory, is created if missing.
 * @throws std::invalid_argument if `config.shot_id` has a forbidden character
 *         or `kind` is not a BackendKind value.
 * @throws std::runtime_error if the back-end is not built in or the output
 *         directory cannot be created.
 */
std::unique_ptr<IOControllerBase> makeIOController(BackendKind kind, OpenMode mode, const IOConfig& config);

}  // namespace funtides::io

#endif  // FUNTIDES_IO_INCLUDE_IO_CONTROLLER_BASE_H_
