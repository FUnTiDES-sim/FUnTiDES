#ifndef FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_VIEW_H_
#define FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_VIEW_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "data_type.h"
#include "wavefield_view.h"
#include "wavefield_view_backward_acoustic.h"
#include "wavefield_view_backward_elastic.h"
#include "wavefield_view_forward_acoustic.h"
#include "wavefield_view_forward_elastic.h"

namespace py = pybind11;

namespace gradient {

void bind_wavefield_view_base(py::module_& m) {
  // Bind WavefieldView (base class)
  py::class_<WavefieldView, std::shared_ptr<WavefieldView>>(m, "WavefieldView").def("print", &WavefieldView::print);
}

void bind_wavefield_view_forward_acoustic(py::module_& m) {
  // Bind WavefieldViewForwardAcoustic (inherits from WavefieldView)
  py::class_<WavefieldViewForwardAcoustic, WavefieldView, std::shared_ptr<WavefieldViewForwardAcoustic>>(
      m, "WavefieldViewForwardAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(), py::arg("pn"))
      .def("print", &WavefieldViewForwardAcoustic::print);
}

void bind_wavefield_view_backward_acoustic(py::module_& m) {
  // Bind WavefieldViewBackwardAcoustic (inherits from WavefieldView)
  py::class_<WavefieldViewBackwardAcoustic, WavefieldView, std::shared_ptr<WavefieldViewBackwardAcoustic>>(
      m, "WavefieldViewBackwardAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("qn"), py::arg("qn_prev"), py::arg("qn_prev_prev"))
      .def("print", &WavefieldViewBackwardAcoustic::print);
}

void bind_wavefield_view_forward_elastic(py::module_& m) {
  // Bind WavefieldViewForwardElastic (inherits from WavefieldView)
  py::class_<WavefieldViewForwardElastic, WavefieldView, std::shared_ptr<WavefieldViewForwardElastic>>(
      m, "WavefieldViewForwardElastic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("ux_n"), py::arg("uy_n"), py::arg("uz_n"))
      .def("print", &WavefieldViewForwardElastic::print);
}

void bind_wavefield_view_backward_elastic(py::module_& m) {
  // Bind WavefieldViewBackwardElastic (inherits from WavefieldView)
  py::class_<WavefieldViewBackwardElastic, WavefieldView, std::shared_ptr<WavefieldViewBackwardElastic>>(
      m, "WavefieldViewBackwardElastic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("ux_n"), py::arg("uy_n"), py::arg("uz_n"), py::arg("ux_dt2"), py::arg("uy_dt2"), py::arg("uz_dt2"))
      .def("print", &WavefieldViewBackwardElastic::print);
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_VIEW_H_