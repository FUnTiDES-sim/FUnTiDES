#ifndef FUNTIDES_SOLVER_FE_DG_PYWRAP_INCLUDE_BINDINGS_DG_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_DG_PYWRAP_INCLUDE_BINDINGS_DG_WAVEFIELD_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "dg_wavefield_acoustic.h"

namespace py = pybind11;

namespace solver {
namespace fe {

void bind_dg_wavefield_acoustic(py::module_ &m) {
  py::class_<DGWavefieldAcoustic, std::shared_ptr<DGWavefieldAcoustic>>(m, "DGWavefieldAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("pn_prev"), py::arg("pn_curr"))
      .def("swap", &DGWavefieldAcoustic::swap)
      .def(
          "get_current_field",
          [](const DGWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const DGWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getPreviousField(i));
          },
          py::arg("i"))
      .def("print", &DGWavefieldAcoustic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PYWRAP_INCLUDE_BINDINGS_DG_WAVEFIELD_H_
