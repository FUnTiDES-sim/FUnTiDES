#ifndef FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_
#define FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <string>

#include "common_macros.h"
#include "data_type.h"
#include "rhs.h"
#include "rhs_acoustic.h"
#include "rhs_elastic.h"

namespace py = pybind11;

namespace solver
{
namespace fe
{

void bind_rhs_base(py::module_ &m)
{
  // Bind Rhs (base class)
  py::class_<Rhs, std::shared_ptr<Rhs>>(m, "Rhs").def("print", &Rhs::print);
}

void bind_rhs_acoustic(py::module_ &m)
{
  // Automatically deduce the underlying scalar types from your macros
  using real_type = typename ARRAY_REAL_VIEW::value_type;
  using int_type = typename VECTOR_INT_VIEW::value_type;

  // Bind RhsAcoustic (inherits from Rhs)
  py::class_<RhsAcoustic, Rhs, std::shared_ptr<RhsAcoustic>>(m, "RhsAcoustic")
      .def(
          py::init([](py::array_t<real_type> term_py,
                      py::array_t<int_type> element_py,
                      py::array_t<real_type> weights_py) {

              // Helper lambda for 2D Real arrays
              auto wrap_and_copy_2d_real = [](py::array_t<real_type> arr, const std::string& name) {
                  auto buf = arr.request();
                  Kokkos::View<real_type**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                      h_view((real_type*)buf.ptr, buf.shape[0], buf.shape[1]);
                  ARRAY_REAL_VIEW d_view(name, buf.shape[0], buf.shape[1]);
                  Kokkos::deep_copy(d_view, h_view);
                  return d_view;
              };

              // Helper lambda for 1D Int arrays
              auto wrap_and_copy_1d_int = [](py::array_t<int_type> arr, const std::string& name) {
                  auto buf = arr.request();
                  Kokkos::View<int_type*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                      h_view((int_type*)buf.ptr, buf.shape[0]);
                  VECTOR_INT_VIEW d_view(name, buf.shape[0]);
                  Kokkos::deep_copy(d_view, h_view);
                  return d_view;
              };

              return new RhsAcoustic(
                  wrap_and_copy_2d_real(term_py, "rhs_term"),
                  wrap_and_copy_1d_int(element_py, "rhs_element"),
                  wrap_and_copy_2d_real(weights_py, "rhs_weights")
              );
          }),
          py::arg("term"), py::arg("element"), py::arg("weights"))
      .def("print", &RhsAcoustic::print);
}

void bind_rhs_elastic(py::module_ &m)
{
  using real_type = typename ARRAY_REAL_VIEW::value_type;
  using int_type = typename VECTOR_INT_VIEW::value_type;

  // Bind RhsElastic (inherits from Rhs)
  py::class_<RhsElastic, Rhs, std::shared_ptr<RhsElastic>>(m, "RhsElastic")
      .def(
          py::init([](py::array_t<real_type> termx_py,
                      py::array_t<real_type> termy_py,
                      py::array_t<real_type> termz_py,
                      py::array_t<int_type> element_py,
                      py::array_t<real_type> weights_py) {

              auto wrap_and_copy_2d_real = [](py::array_t<real_type> arr, const std::string& name) {
                  auto buf = arr.request();
                  Kokkos::View<real_type**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                      h_view((real_type*)buf.ptr, buf.shape[0], buf.shape[1]);
                  ARRAY_REAL_VIEW d_view(name, buf.shape[0], buf.shape[1]);
                  Kokkos::deep_copy(d_view, h_view);
                  return d_view;
              };

              auto wrap_and_copy_1d_int = [](py::array_t<int_type> arr, const std::string& name) {
                  auto buf = arr.request();
                  Kokkos::View<int_type*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                      h_view((int_type*)buf.ptr, buf.shape[0]);
                  VECTOR_INT_VIEW d_view(name, buf.shape[0]);
                  Kokkos::deep_copy(d_view, h_view);
                  return d_view;
              };

              return new RhsElastic(
                  wrap_and_copy_2d_real(termx_py, "rhs_term_x"),
                  wrap_and_copy_2d_real(termy_py, "rhs_term_y"),
                  wrap_and_copy_2d_real(termz_py, "rhs_term_z"),
                  wrap_and_copy_1d_int(element_py, "rhs_element"),
                  wrap_and_copy_2d_real(weights_py, "rhs_weights")
              );
          }),
          py::arg("termx"), py::arg("termy"), py::arg("termz"),
          py::arg("element"), py::arg("weights"))
      .def("print", &RhsElastic::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_RHS_H_
