#ifndef FUNTIDES_SOLVER_FE_DG_SEM_PYWRAP_INCLUDE_BINDINGS_DGSEM_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_PYWRAP_INCLUDE_BINDINGS_DGSEM_WAVEFIELD_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "dg-sem_wavefield_acoustic.h"

namespace py = pybind11;

namespace solver {
namespace fe {

void bind_dgsem_wavefield_acoustic(py::module_ &m) {
  py::class_<DGSEMWavefieldAcoustic, std::shared_ptr<DGSEMWavefieldAcoustic>>(m, "DGSEMWavefieldAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>>(),
           py::arg("pn_dg_prev"), py::arg("pn_dg_curr"), py::arg("pn_sem_prev"), py::arg("pn_sem_curr"))
      .def("swap", &DGSEMWavefieldAcoustic::swap)
      .def(
          "get_dg_current_field",
          [](const DGSEMWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getDGCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_sem_current_field",
          [](const DGSEMWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getSEMCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_dg_previous_field",
          [](const DGSEMWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<arrayReal>(self.getDGPreviousField(i));
          },
          py::arg("i"))
      .def(
          "get_sem_previous_field",
          [](const DGSEMWavefieldAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getSEMPreviousField(i));
          },
          py::arg("i"))
      .def("print", &DGSEMWavefieldAcoustic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_SEM_PYWRAP_INCLUDE_BINDINGS_DGSEM_WAVEFIELD_H_
