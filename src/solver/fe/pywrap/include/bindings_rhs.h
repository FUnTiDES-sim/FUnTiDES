#ifndef FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_
#define FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "data_type.h"
#include "rhs.h"
#include "rhs_acoustic.h"
#include "rhs_acoustoelastic.h"
#include "rhs_elastic.h"

namespace py = pybind11;

namespace solver {
namespace fe {

void bind_rhs_base(py::module_& m) {
  // Bind Rhs (base class)
  py::class_<Rhs, std::shared_ptr<Rhs>>(m, "Rhs").def("print", &Rhs::print);
}

void bind_rhs_acoustic(py::module_& m) {
  // Bind RhsAcoustic (inherits from Rhs)
  py::class_<RhsAcoustic, Rhs, std::shared_ptr<RhsAcoustic>>(m, "RhsAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<vectorInt>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("term"), py::arg("element"), py::arg("weights"))
      .def("print", &RhsAcoustic::print);
}

void bind_rhs_elastic(py::module_& m) {
  // Bind RhsElastic (inherits from Rhs)
  py::class_<RhsElastic, Rhs, std::shared_ptr<RhsElastic>>(m, "RhsElastic")
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<arrayReal>, Kokkos::Experimental::python_view_type_t<arrayReal>,
               Kokkos::Experimental::python_view_type_t<arrayReal>, Kokkos::Experimental::python_view_type_t<vectorInt>,
               Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("termx"), py::arg("termy"), py::arg("termz"), py::arg("element"), py::arg("weights"))
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<arrayReal>, Kokkos::Experimental::python_view_type_t<arrayReal>,
               Kokkos::Experimental::python_view_type_t<arrayReal>, Kokkos::Experimental::python_view_type_t<vectorInt>,
               Kokkos::Experimental::python_view_type_t<arrayReal>, Kokkos::Experimental::python_view_type_t<arrayReal>,
               Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("termx"), py::arg("termy"), py::arg("termz"), py::arg("element"), py::arg("weightsx"),
           py::arg("weightsy"), py::arg("weightsz"))
      .def("print", &RhsElastic::print);
}

void bind_rhs_acoustoelastic(py::module_& m) {
  py::class_<RhsAcoustoElastic, Rhs, std::shared_ptr<RhsAcoustoElastic>>(m, "RhsAcoustoElastic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,  // acoustic_term
                    Kokkos::Experimental::python_view_type_t<vectorInt>,  // element
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // weights
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_termx
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_termy
                    Kokkos::Experimental::python_view_type_t<arrayReal>   // elastic_termz
                    >(),
           py::arg("acoustic_term"), py::arg("element"), py::arg("weights"), py::arg("elastic_termx"),
           py::arg("elastic_termy"), py::arg("elastic_termz"))
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,  // acoustic_term
                    Kokkos::Experimental::python_view_type_t<vectorInt>,  // element
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // weights
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_termx
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_termy
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_termz
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_weightsx
                    Kokkos::Experimental::python_view_type_t<arrayReal>,  // elastic_weightsy
                    Kokkos::Experimental::python_view_type_t<arrayReal>   // elastic_weightsz
                    >(),
           py::arg("acoustic_term"), py::arg("element"), py::arg("weights"), py::arg("elastic_termx"),
           py::arg("elastic_termy"), py::arg("elastic_termz"), py::arg("elastic_weightsx"), py::arg("elastic_weightsy"),
           py::arg("elastic_weightsz"))
      .def("print", &RhsAcoustoElastic::print);
}
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_