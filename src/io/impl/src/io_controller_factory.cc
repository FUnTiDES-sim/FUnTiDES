#include <stdexcept>

#include "io_controller_base.h"
#include "posix_io_controller.h"

namespace funtides::io {

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
