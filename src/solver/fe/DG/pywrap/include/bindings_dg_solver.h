#ifndef FUNTIDES_SOLVER_FE_DG_PYWRAP_INCLUDE_BINDINGS_DG_SOLVER_H_
#define FUNTIDES_SOLVER_FE_DG_PYWRAP_INCLUDE_BINDINGS_DG_SOLVER_H_
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "dg_solver_data.h"
#include "rhs_acoustic.h"

namespace py = pybind11;

namespace solver {
namespace fe {

void bind_dg_acoustic_data(py::module_ &m) {
  py::class_<DGsolverDataAcoustic, Solver::DataStruct, std::shared_ptr<DGsolverDataAcoustic>>(m, "DGsolverDataAcoustic")
      .def(py::init<const DGWavefieldAcoustic &, const RhsAcoustic &>(), py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &DGsolverDataAcoustic::swapWavefields)
      .def("print", &DGsolverDataAcoustic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PYWRAP_INCLUDE_BINDINGS_DG_SOLVER_H_
