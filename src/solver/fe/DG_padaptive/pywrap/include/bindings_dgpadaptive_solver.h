#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_PYWRAP_INCLUDE_BINDINGS_DGPADAPTIVE_SOLVER_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_PYWRAP_INCLUDE_BINDINGS_DGPADAPTIVE_SOLVER_H_
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "dg_padaptive_rhs_acoustic.h"
#include "dg_padaptive_solver_data.h"
#include "rhs.h"

namespace py = pybind11;

namespace solver {
namespace fe {

/**
 * @brief Registers the p-adaptive acoustic source term in a Python module.
 *
 * The Python constructor takes, in order, the pMin and pMax acoustic terms, the source element
 * indices, and the pMin and pMax weights. Exposes `print`.
 *
 * @param[in,out] m Module receiving the class.
 * @todo VERIFY: shapes and layouts of the term and weight arrays expected by the Python constructor.
 */
void bind_dgpadaptive_rhs_acoustic(py::module_ &m) {
  py::class_<DGPAdaptiveRhsAcoustic, Rhs, std::shared_ptr<DGPAdaptiveRhsAcoustic>>(m, "DGPAdaptiveRhsAcoustic")
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<arrayReal>, Kokkos::Experimental::python_view_type_t<arrayReal>,
               Kokkos::Experimental::python_view_type_t<vectorInt>, Kokkos::Experimental::python_view_type_t<arrayReal>,
               Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("pmin_acoustic_term"), py::arg("pmax_acoustic_term"), py::arg("element"), py::arg("pmin_weights"),
           py::arg("pmax_weights"))
      .def("print", &DGPAdaptiveRhsAcoustic::print);
}

/**
 * @brief Registers the p-adaptive solver data (wavefield and source) in a Python module.
 *
 * The Python constructor takes a wavefield and a source term. Exposes `swap_wavefields` and `print`.
 *
 * @param[in,out] m Module receiving the class.
 */
void bind_dgpadaptive_acoustic_data(py::module_ &m) {
  py::class_<DGPAdaptiveSolverData, Solver::DataStruct, std::shared_ptr<DGPAdaptiveSolverData>>(m,
                                                                                                "DGPAdaptiveSolverData")
      .def(py::init<const DGPAdaptiveWavefieldAcoustic &, const DGPAdaptiveRhsAcoustic &>(), py::arg("wavefield"),
           py::arg("rhs"))
      .def("swap_wavefields", &DGPAdaptiveSolverData::swapWavefields)
      .def("print", &DGPAdaptiveSolverData::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_PYWRAP_INCLUDE_BINDINGS_DGPADAPTIVE_SOLVER_H_
