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

namespace solver
{
namespace fe
{

void bind_rhs_base(py::module_& m)
{
  // Bind Rhs (base class)
  py::class_<Rhs, std::shared_ptr<Rhs>>(m, "Rhs").def("print", &Rhs::print);
}

void bind_rhs_acoustic(py::module_& m)
{
  // Bind RhsAcoustic (inherits from Rhs)
  py::class_<RhsAcoustic, Rhs, std::shared_ptr<RhsAcoustic>>(m, "RhsAcoustic")
      .def(
          py::init<Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<VECTOR_INT_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>>(),
          py::arg("term"), py::arg("element"), py::arg("weights"))
      .def("print", &RhsAcoustic::print);
}

void bind_rhs_elastic(py::module_& m)
{
  // Bind RhsElastic (inherits from Rhs)
  py::class_<RhsElastic, Rhs, std::shared_ptr<RhsElastic>>(m, "RhsElastic")
      .def(
          py::init<Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<VECTOR_INT_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>>(),
          py::arg("termx"), py::arg("termy"), py::arg("termz"),
          py::arg("element"), py::arg("weights"))
      .def("print", &RhsElastic::print);
}

void bind_rhs_acoustoelastic(py::module_& m)
{
  using real_type = typename ARRAY_REAL_VIEW::value_type;
  using int_type = typename VECTOR_INT_VIEW::value_type;

  // Bind RhsAcoustoElastic (inherits from Rhs)
  py::class_<RhsAcoustoElastic, Rhs, std::shared_ptr<RhsAcoustoElastic>>(
      m, "RhsAcoustoElastic")
      .def(py::init([](Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>
                           acoustic_term_py,  // CORRIGÉ
                       Kokkos::Experimental::python_view_type_t<VECTOR_INT_VIEW>
                           element_py,  // CORRIGÉ
                       Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>
                           weights_py,  // CORRIGÉ
                       Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>
                           elastic_termx_py,  // CORRIGÉ
                       Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>
                           elastic_termy_py,  // CORRIGÉ
                       Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>
                           elastic_termz_py) {  // CORRIGÉ
             // Les lambdas internes doivent aussi utiliser les types pybind11
             // corrects Souvent on utilise simplement py::array_t<real_type>
             // pour plus de clarté

             auto wrap_and_copy_2d_real = [](py::array_t<real_type> arr,
                                             const std::string& name) {
               auto buf = arr.request();
               if (buf.ndim != 2)
                 throw std::runtime_error(name + " must be 2D");

               Kokkos::View<real_type**, Kokkos::LayoutRight, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                   h_view((real_type*)buf.ptr, buf.shape[0], buf.shape[1]);

               ARRAY_REAL_VIEW d_view(name, buf.shape[0], buf.shape[1]);
               Kokkos::deep_copy(d_view, h_view);
               return d_view;
             };

             auto wrap_and_copy_1d_int = [](py::array_t<int_type> arr,
                                            const std::string& name) {
               auto buf = arr.request();
               if (buf.ndim != 1)
                 throw std::runtime_error(name + " must be 1D");

               Kokkos::View<int_type*, Kokkos::LayoutRight, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                   h_view((int_type*)buf.ptr, buf.shape[0]);

               VECTOR_INT_VIEW d_view(name, buf.shape[0]);
               Kokkos::deep_copy(d_view, h_view);
               return d_view;
             };

             return new RhsAcoustoElastic(
                 wrap_and_copy_2d_real(acoustic_term_py, "rhs_acoustic_term"),
                 wrap_and_copy_1d_int(element_py, "rhs_element"),
                 wrap_and_copy_2d_real(weights_py, "rhs_weights"),
                 wrap_and_copy_2d_real(elastic_termx_py, "rhs_elastic_term_x"),
                 wrap_and_copy_2d_real(elastic_termy_py, "rhs_elastic_term_y"),
                 wrap_and_copy_2d_real(elastic_termz_py, "rhs_elastic_term_z"));
           }),
           py::arg("acoustic_term"), py::arg("element"), py::arg("weights"),
           py::arg("elastic_termx"), py::arg("elastic_termy"),
           py::arg("elastic_termz"))
      .def("print", &RhsAcoustoElastic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_