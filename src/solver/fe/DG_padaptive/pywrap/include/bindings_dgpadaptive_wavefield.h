#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_PYWRAP_INCLUDE_BINDINGS_DGPADAPTIVE_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_PYWRAP_INCLUDE_BINDINGS_DGPADAPTIVE_WAVEFIELD_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "dg_padaptive_wavefield_acoustic.h"

namespace py = pybind11;

namespace solver {
namespace fe {

void bind_dgpadaptive_wavefield_acoustic(py::module_ &m) {
  py::class_<DGPAdaptiveWavefieldAcoustic, std::shared_ptr<DGPAdaptiveWavefieldAcoustic>>(
      m, "DGPAdaptiveWavefieldAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("pn_pmin_prev"), py::arg("pn_pmin_curr"), py::arg("pn_pmax_prev"), py::arg("pn_pmax_curr"))
      .def("swap", &DGPAdaptiveWavefieldAcoustic::swap)
      .def(
          "get_pmin_current_field",
          [](const DGPAdaptiveWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getPMinCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_pmax_current_field",
          [](const DGPAdaptiveWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getPMaxCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_pmin_previous_field",
          [](const DGPAdaptiveWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getPMinPreviousField(i));
          },
          py::arg("i"))
      .def(
          "get_pmax_previous_field",
          [](const DGPAdaptiveWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getPMaxPreviousField(i));
          },
          py::arg("i"))
      .def("print", &DGPAdaptiveWavefieldAcoustic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_PYWRAP_INCLUDE_BINDINGS_DGPADAPTIVE_WAVEFIELD_H_
