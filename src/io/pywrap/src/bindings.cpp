#include <pybind11/pybind11.h>

#include "bindings_io.h"

namespace py = pybind11;

/**
 * @brief Python module exposing snapshot I/O.
 *
 * The I/O API is not templated: snapshots use the solver's real type, so a
 * single set of classes is registered. The ADIOS2 backend is always listed in
 * BackendKind; is_backend_available() tells whether this build has it.
 */
PYBIND11_MODULE(io, m) {
  m.attr("__name__") = "pyfuntides.io";
  m.doc() = "Wavefield snapshot I/O (POSIX and optional ADIOS2 backends).";

  io_bindings::bind_io_enums(m);
  io_bindings::bind_io_config(m);
  io_bindings::bind_io_controller(m);
}
