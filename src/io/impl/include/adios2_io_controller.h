#ifndef FUNTIDES_IO_INCLUDE_ADIOS2_IO_CONTROLLER_H_
#define FUNTIDES_IO_INCLUDE_ADIOS2_IO_CONTROLLER_H_

#include <memory>

#include "io_controller_base.h"

namespace funtides::io {

/**
 * @brief I/O backend storing the snapshots of one rank as the steps of one
 * ADIOS2 BP5 file, with the file accesses done by a background host thread.
 *
 * File layout: <output_dir>[/<shot_id>]/<prefix>_r<rank>.bp
 * Step `i` of the file holds snapshot `i`, in a 1D variable named "field".
 * `local_dims`, `global_dims`, `start_offsets` and the rank are stored as
 * attributes, so a post-processing tool can reassemble the global grid.
 *
 * ASYNCHRONY. Each controller owns one I/O thread, and that thread is the only
 * one that calls ADIOS2.
 * - Write: writeSnapshot() copies `field` into an internal host staging
 *   buffer, queues it and returns; the thread writes it. The caller's buffer
 *   is free as soon as writeSnapshot() returns (stronger than the base
 *   contract). With `async_snapshots == false`, writeSnapshot() waits for the
 *   write to land before returning. When every staging buffer is in flight,
 *   writeSnapshot() blocks until one is free (back-pressure: memory use is
 *   bounded by kWriteBuffers snapshots).
 * - Read: after readSnapshot(i), the thread prefetches the next snapshots in
 *   the direction of the access pattern (i+1, i+2... or i-1, i-2... for a
 *   backward sweep, e.g. the adjoint pass of an RTM). A readSnapshot() that
 *   hits a prefetched snapshot only costs a host memcpy.
 *
 * ERRORS. An error raised by the I/O thread is rethrown by the next call to
 * writeSnapshot(), readSnapshot(), flush() or close(); the controller stays
 * failed afterwards.
 *
 * MPI. Every rank writes its own file with the serial ADIOS2 API: the I/O
 * thread never makes an MPI call, so no MPI_THREAD_MULTIPLE is required and
 * close() is not collective in practice.
 *
 * Every snapshot of a controller must have the same size.
 */
class Adios2IOController final : public IOControllerBase {
 public:
  /// Number of host staging buffers of a write controller.
  static constexpr std::size_t kWriteBuffers = 2;
  /// Number of snapshots a read controller loads ahead of the caller.
  static constexpr std::size_t kPrefetchDepth = 2;

  /**
   * @brief Opens the BP file and starts the I/O thread.
   *
   * Returns once the file is open, so open errors are thrown here.
   * @throws std::invalid_argument if `config.shot_id` has a forbidden character.
   * @throws std::runtime_error if the directory cannot be created (write), the
   *         file does not exist or holds no "field" variable of the right
   *         scalar type (read), or ADIOS2 fails to open it.
   */
  Adios2IOController(OpenMode mode, const IOConfig& config);

  /// @brief Closes the controller; errors are reported on stderr only.
  ~Adios2IOController() override;

  Adios2IOController(const Adios2IOController&) = delete;
  Adios2IOController& operator=(const Adios2IOController&) = delete;

  /**
   * @brief Copies `field` into a staging buffer and queues it for writing.
   * @throws std::runtime_error if the size differs from the first snapshot,
   *         `field` is not contiguous, or a previous write failed.
   */
  void writeSnapshot(const HostVectorReal& field) override;

  /**
   * @brief Copies snapshot `index` into `field`, loading it first unless it
   * was prefetched, then starts prefetching the next ones.
   * @throws std::runtime_error if `index` is out of range or the size of
   *         `field` differs from the stored one.
   */
  void readSnapshot(const HostVectorReal& field, std::size_t index) override;

  /// @brief Write mode: waits until every queued snapshot is written. Read mode: no-op.
  void flush() override;

  /// @brief Drains the write queue, stops the I/O thread and closes the file. Idempotent.
  void close() override;

 private:
  class Impl;
  void requireMode(OpenMode expected, const char* what) const;

  OpenMode mode_;
  bool closed_{false};
  std::unique_ptr<Impl> impl_;
};

}  // namespace funtides::io

#endif  // FUNTIDES_IO_INCLUDE_ADIOS2_IO_CONTROLLER_H_
