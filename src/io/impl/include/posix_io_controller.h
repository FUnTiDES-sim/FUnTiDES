#ifndef FUNTIDES_IO_INCLUDE_POSIX_IO_CONTROLLER_H_
#define FUNTIDES_IO_INCLUDE_POSIX_IO_CONTROLLER_H_

#include <cstddef>
#include <string>

#include "io_controller_base.h"

namespace funtides::io {

/**
 * @brief I/O backend writing one plain binary file per snapshot and per rank.
 *
 * File layout: <output_dir>[/<shot_id>]/<prefix>_snap_<index>_r<rank>.bin
 * Each file starts with a fixed-size self-describing header, so a reader can
 * validate extents without external metadata. Writes are synchronous: the
 * caller's buffers are free as soon as the call returns, and `async_snapshots`
 * is ignored.
 *
 * MPI: ranks write independent files; nothing is merged. The global shape and
 * this rank's offset are recorded in each header so a post-processing tool can
 * reassemble.
 */
class PosixIOController final : public IOControllerBase {
 public:
  /**
   * @brief Opens the controller in the given mode.
   * @param[in] mode Whether snapshots are written or read.
   * @param[in] config Output directory, file prefix, shot id and layout settings.
   */
  PosixIOController(OpenMode mode, const IOConfig& config);

  /// @brief Releases the controller.
  ~PosixIOController() override;

  /**
   * @brief Writes the next snapshot of this rank to a new file.
   * @param[in] field Host data of the local subdomain.
   */
  void writeSnapshot(const HostVectorReal& field) override;

  /**
   * @brief Reads the snapshot with the given index into a host vector.
   * @param[in] field Destination buffer; the header extents are validated against it.
   * @param[in] index Snapshot index, as used in the file name.
   * @todo VERIFY: is `field` really written to (const reference) and what happens on extent mismatch?
   */
  void readSnapshot(const HostVectorReal& field, std::size_t index) override;

  /// @brief Flushes pending writes to disk.
  void flush() override;

  /// @brief Closes the current file; the controller is unusable afterwards.
  void close() override;

 private:
  /// Directory holding this shot's snapshots.
  std::string outputDir() const;
  /// Full path of the file of this rank for snapshot `index`.
  std::string snapshotPath(std::size_t index) const;
  /// Throws if the controller was not opened in mode `expected`; `what` names the caller in the message.
  void requireMode(OpenMode expected, const char* what) const;

  OpenMode mode_;              ///< Mode chosen at construction.
  std::size_t next_index_{0};  ///< Index given to the next written snapshot.
  bool closed_{false};         ///< True once close() has been called.
};

}  // namespace funtides::io

#endif  // FUNTIDES_IO_INCLUDE_POSIX_IO_CONTROLLER_H_
