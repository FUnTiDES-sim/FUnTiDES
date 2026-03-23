#ifndef SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_GRADIENTS_H_
#define SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_GRADIENTS_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "gradients.h"
#include "gradients_acoustic.h"

namespace py = pybind11;

namespace solver
{
namespace fe
{

void bind_gradients_base(py::module_ &m)
{
  // Bind Gradients (base class)
  py::class_<Gradients, std::shared_ptr<Gradients>>(m, "Gradients")
      .def("print", &Gradients::print);
}

void bind_gradients_acoustic(py::module_ &m)
{
  // Bind GradientsAcoustic (inherits from Gradients)
  py::class_<GradientsAcoustic, Gradients, std::shared_ptr<GradientsAcoustic>>(
      m, "GradientsAcoustic")
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("grad_kappa"), py::arg("grad_buoyancy"))
      .def("print", &GradientsAcoustic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_GRADIENTS_H_
