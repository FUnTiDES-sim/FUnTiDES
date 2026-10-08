#include <stdexcept>

#include "io_controller_base.h"
#include "posix_io_controller.h"

#ifdef FUNTIDES_IO_HAVE_ADIOS2
#include "adios2_io_controller.h"
#endif

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
#ifdef FUNTIDES_IO_HAVE_ADIOS2
      return std::make_unique<Adios2IOController>(mode, config);
#else
      throw std::runtime_error(
          "funtides::io: ADIOS2 backend not built in (configure with -DFUNTIDES_ENABLE_ADIOS2=ON)");
#endif
  }
  throw std::invalid_argument("funtides::io: unknown backend");
}

bool isBackendAvailable(BackendKind kind) noexcept {
  switch (kind) {
    case BackendKind::kPosix:
      return true;
    case BackendKind::kAdios2:
#ifdef FUNTIDES_IO_HAVE_ADIOS2
      return true;
#else
      return false;
#endif
  }
  return false;
}

}  // namespace funtides::io
