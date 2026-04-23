#ifndef FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_GRADIENT_H_
#define FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_GRADIENT_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "data_type.h"
#include "gradient.h"
#include "gradient_acoustic.h"
#include "gradient_elastic.h"

namespace py = pybind11;

namespace gradient {

void bind_gradient_base(py::module_& m) {
  // Bind Gradient (base class)
  py::class_<Gradient, std::shared_ptr<Gradient>>(m, "Gradient").def("print", &Gradient::print);
}

void bind_gradient_acoustic(py::module_& m) {
  // Bind GradientAcoustic (inherits from Gradient)
  py::class_<GradientAcoustic, Gradient, std::shared_ptr<GradientAcoustic>>(m, "GradientAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("grad_kappa"), py::arg("grad_buoyancy"))
      .def("print", &GradientAcoustic::print);
}

void bind_gradient_elastic(py::module_& m) {
  // Bind GradientElastic (inherits from Gradient)
  py::class_<GradientElastic, Gradient, std::shared_ptr<GradientElastic>>(m, "GradientElastic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
                    Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("grad_rho"), py::arg("grad_lambda"), py::arg("grad_mu"))
      .def("print", &GradientElastic::print);
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_GRADIENT_H_