#ifndef FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "wavefield.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace py = pybind11;

namespace solver
{
namespace fe
{

void bind_wavefield_base(py::module_& m)
{
  // Bind Wavefield (base class)
  py::class_<Wavefield, std::shared_ptr<Wavefield>>(m, "Wavefield")
      .def("swap", &Wavefield::swap)
      .def(
          "swap_with_rotation",
          [](Wavefield& self,
             Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>&
                 prevPrevBuffer) {
            VECTOR_REAL_VIEW view = prevPrevBuffer;
            self.swapWithRotation(view, 0);  // field index 0 for pressure
            prevPrevBuffer = view;
          },
          py::arg("prev_prev_buffer"))
      .def(
          "get_current_field",
          [](const Wavefield& self, int i) {
            return Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>(
                self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const Wavefield& self, int i) {
            return Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>(
                self.getPreviousField(i));
          },
          py::arg("i"))
      .def("print", &Wavefield::print);
}

void bind_wavefield_acoustic(py::module_& m)
{
  // Bind WavefieldAcoustic (inherits from Wavefield)
  py::class_<WavefieldAcoustic, Wavefield, std::shared_ptr<WavefieldAcoustic>>(
      m, "WavefieldAcoustic")
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("pn_global_prev"), py::arg("pn_global_curr"))
      .def("swap", &WavefieldAcoustic::swap)
      .def(
          "swap_with_rotation",
          [](WavefieldAcoustic& self,
             Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>&
                 prevPrevBuffer) {
            VECTOR_REAL_VIEW view = prevPrevBuffer;
            self.swapWithRotation(view, 0);  // field index 0 for pressure
            prevPrevBuffer = view;
          },
          py::arg("prev_prev_buffer"))
      .def(
          "get_current_field",
          [](const WavefieldAcoustic& self, int i) {
            return Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>(
                self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const WavefieldAcoustic& self, int i) {
            return Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>(
                self.getPreviousField(i));
          },
          py::arg("i"))
      .def("print", &WavefieldAcoustic::print);
}

void bind_wavefield_elastic(py::module_& m)
{
  // Bind WavefieldElastic (inherits from Wavefield)
  py::class_<WavefieldElastic, Wavefield, std::shared_ptr<WavefieldElastic>>(
      m, "WavefieldElastic")
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("uxn_global_prev"), py::arg("uxn_global_curr"),
           py::arg("uyn_global_prev"), py::arg("uyn_global_curr"),
           py::arg("uzn_global_prev"), py::arg("uzn_global_curr"))
      .def("swap", &WavefieldElastic::swap)
      .def(
          "swap_with_rotation",
          [](WavefieldElastic& self,
             Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>& uxPP,
             Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>& uyPP,
             Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>& uzPP) {
            VECTOR_REAL_VIEW vux = uxPP, vuy = uyPP, vuz = uzPP;
            self.swapWithRotation(vux, 0);  // ux component
            self.swapWithRotation(vuy, 1);  // uy component
            self.swapWithRotation(vuz, 2);  // uz component
            uxPP = vux;
            uyPP = vuy;
            uzPP = vuz;
          },
          py::arg("ux_prev_prev"), py::arg("uy_prev_prev"),
          py::arg("uz_prev_prev"))
      .def(
          "get_current_field",
          [](const WavefieldElastic& self, int i) {
            return Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>(
                self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const WavefieldElastic& self, int i) {
            return Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>(
                self.getPreviousField(i));
          },
          py::arg("i"))
      .def("print", &WavefieldElastic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_H_