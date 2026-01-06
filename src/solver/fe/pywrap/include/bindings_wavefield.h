#pragma once

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

void bind_wavefield_base(py::module_ &m)
{
  // Bind Wavefield (base class)
  py::class_<Wavefield, std::shared_ptr<Wavefield>>(m, "Wavefield")
      .def("swap", &Wavefield::swap)
      .def("print", &Wavefield::print);
}

void bind_wavefield_acoustic(py::module_ &m)
{
  // Bind WavefieldAcoustic (inherits from Wavefield)
  py::class_<WavefieldAcoustic, Wavefield, std::shared_ptr<WavefieldAcoustic>>(
      m, "WavefieldAcoustic")
      .def(py::init<
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("pn_global_prev"), py::arg("pn_global_curr"))
      .def("swap", &WavefieldAcoustic::swap)
      .def("print", &WavefieldAcoustic::print);
}

void bind_wavefield_elastic(py::module_ &m)
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
      .def("print", &WavefieldElastic::print);
}

}  // namespace fe
}  // namespace solver
