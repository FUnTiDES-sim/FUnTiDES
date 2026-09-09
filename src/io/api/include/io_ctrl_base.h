#ifndef FUNTIDES_IO_INCLUDE_IO_CONTROLLER_BASE_H_
#define FUNTIDES_IO_INCLUDE_IO_CONTROLLER_BASE_H_

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "data_type.h"
#include "distributed_ctx.h"

namespace funtides::io {

using HostVectorReal = vectorReal::host_mirror_type;
using HostArrayReal = arrayReal::host_mirror_type;

struct IOConfig {
  utils::DistributedContext ctx{};
  std::string output_dir{"."};
  std::string prefix{"funtides"};

  std::vector<std::size_t> global_dims;    ///< Global snapshot grid.
  std::vector<std::size_t> start_offsets;  ///< This rank's offset in it.
  std::vector<std::size_t> local_dims;     ///< This rank's chunk.

  std::size_t nb_iter{0};
  std::size_t nb_receiver{0};

  bool async_snapshots{false};
};

/**
 * @brief Backend-agnostic I/O interface for SEM simulation output.
 *
 * HOST ONLY. Never capture an instance inside a Kokkos kernel.
 *
 * THREADING: not thread-safe. A single thread drives the controller.
 *
 * FILE LAYOUT: multi-dimensional arrays are stored row-major regardless of the
 * build's `Layout`, so files written by a CPU build and a GPU build are
 * interchangeable.
 */
class IOControllerBase {
 public:
  virtual ~IOControllerBase() = default;

  /**
   * @brief Writes one snapshot of the pressure/displacement field.
   *
   * BUFFER CONTRACT: `field` must stay valid AND unmodified until the next
   * call to writeSnapshot(), flush(), or close(). Asynchronous backends only
   * read the memory later. Use two host mirrors alternately if the solver
   * needs to keep producing.
   */
  virtual void writeSnapshot(const HostVectorReal& field, int timestep, float time) = 0;

  /**
   * @brief Writes receiver traces and their coordinates.
   *
   * @param traces  Shape {nb_receiver, nb_iter}.
   * @param coords  Shape {nb_receiver, 3}.
   *
   * Synchronous: both buffers are free as soon as this returns.
   */
  virtual void writeReceivers(const HostArrayReal& traces, const HostArrayReal& coords) = 0;

  /**
   * @brief Reads snapshot number `index` into a caller-allocated view.
   * @throws std::runtime_error if `index` is out of range or the extents of
   *         `field` do not match what is stored.
   */
  virtual void readSnapshot(const HostVectorReal& field, std::size_t index) = 0;

  /// Same contract: `traces` must already be sized {nb_receiver, nb_iter}.
  virtual void readReceivers(const HostArrayReal& traces) = 0;

  /// Blocks until every pending write has landed. Buffers handed to
  /// writeSnapshot() are free afterwards.
  virtual void flush() = 0;

  /// Idempotent. Flushes, then closes the files. Propagates errors, so call it
  /// explicitly: the destructor can only swallow them.
  /// COLLECTIVE — every rank must call it, in the same order.
  virtual void close() = 0;

 protected:
  explicit IOControllerBase(const IOConfig& config) : config_(config) {}

  const IOConfig& config() const { return config_; }

 private:
  IOConfig config_;
};

// ---------------------------------------------------------------------------

enum class BackendKind { kPosix, kAdios2 };
enum class OpenMode { kWrite, kRead };

std::unique_ptr<IOControllerBase> makeIOController(BackendKind kind, OpenMode mode, const IOConfig& config);

}  // namespace funtides::io

#endif  // FUNTIDES_IO_INCLUDE_IO_CONTROLLER_BASE_H_
