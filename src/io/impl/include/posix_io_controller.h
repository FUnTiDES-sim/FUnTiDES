#ifndef FUNTIDES_IO_INCLUDE_POSIX_IO_CONTROLLER_H_
#define FUNTIDES_IO_INCLUDE_POSIX_IO_CONTROLLER_H_

#include <cstddef>
#include <string>
#include <vector>

#include "io_controller_base.h"

namespace funtides::io {

/**
 * @brief Reference I/O backend writing plain binary files through <cstdio>.
 *
 * One file per snapshot and per rank:
 *   <output_dir>/<prefix>_snap_<index>_r<rank>.bin
 *   <output_dir>/<prefix>_receivers.bin
 *
 * Each file starts with a fixed-size self-describing header, so a reader can
 * validate extents without external metadata. Writes are synchronous: the
 * caller's buffers are free as soon as the call returns, and `async_snapshots`
 * is ignored.
 *
 * MPI: ranks write independent files; nothing is merged. The global shape and
 * this rank's offset are recorded in each header so a post-processing tool can
 * reassemble. Receiver traces are replicated data and are written by rank 0
 * only.
 */
class PosixIOController final : public IOControllerBase {
 public:
  PosixIOController(OpenMode mode, const IOConfig& config);
  ~PosixIOController() override;

  void writeSnapshot(const HostVectorReal& field, int timestep, float time) override;
  void writeReceivers(const HostArrayReal& traces, const HostArrayReal& coords) override;

  void readSnapshot(const HostVectorReal& field, std::size_t index) override;
  void readReceivers(const HostArrayReal& traces) override;

  void flush() override;
  void close() override;

 private:
  std::string snapshotPath(std::size_t index) const;
  std::string receiversPath() const;
  void requireMode(OpenMode expected, const char* what) const;

  OpenMode mode_;
  std::size_t next_index_{0};
  bool closed_{false};

  std::vector<float> pack_;
};

}  // namespace funtides::io

#endif  // FUNTIDES_IO_INCLUDE_POSIX_IO_CONTROLLER_H_
