#include <stdexcept>

#include "io_controller_base.h"
#include "posix_io_controller.h"

namespace funtides::io {

/**
 * @brief Creates an I/O controller for the requested back-end.
 * @param[in] kind Back-end to instantiate.
 * @param[in] mode Open mode passed to the controller.
 * @param[in] config I/O configuration passed to the controller.
 * @return Owning pointer to the new controller.
 * @throws std::runtime_error if the requested back-end is not built in.
 * @throws std::invalid_argument if @p kind is not a known back-end.
 */
std::unique_ptr<IOControllerBase> makeIOController(BackendKind kind, OpenMode mode, const IOConfig& config) {
  switch (kind) {
    case BackendKind::kPosix:
      return std::make_unique<PosixIOController>(mode, config);
    case BackendKind::kAdios2:
      throw std::runtime_error("funtides::io: ADIOS2 backend not built in");
  }
  throw std::invalid_argument("funtides::io: unknown backend");
}

}  // namespace funtides::io
