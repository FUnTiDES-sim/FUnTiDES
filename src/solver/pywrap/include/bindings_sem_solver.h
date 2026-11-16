#ifndef FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <memory>

#include "common_macros.h"
#include "rhs_acoustic.h"
#include "rhs_acoustoelastic.h"
#include "rhs_elastic.h"
#include "sem_solver.h"
#include "sem_solver_acoustoelastic.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"
#include "wavefield_acoustoelastic.h"
#include "wavefield_elastic.h"

namespace py = pybind11;

namespace solver {

void bind_acoustic_solver_data(py::module_ &m) {
  py::class_<SEMsolverDataAcoustic, Solver::DataStruct, std::shared_ptr<SEMsolverDataAcoustic>>(m,
                                                                                                "SEMsolverDataAcoustic")
      .def(py::init<const WavefieldAcoustic &, const RhsAcoustic &>(), py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &SEMsolverDataAcoustic::swapWavefields)
      .def(
          "get_current_field",
          [](const SEMsolverDataAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const SEMsolverDataAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPreviousField(i));
          },
          py::arg("i"))
      .def(
          "get_prevprev_field",
          [](const SEMsolverDataAcoustic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPrevPrevField(i));
          },
          py::arg("i"))
      .def("print", &SEMsolverDataAcoustic::print);
}

void bind_elastic_solver_data(py::module_ &m) {
  py::class_<SEMsolverDataElastic, Solver::DataStruct, std::shared_ptr<SEMsolverDataElastic>>(m, "SEMsolverDataElastic")
      .def(py::init<const WavefieldElastic &, const RhsElastic &>(), py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &SEMsolverDataElastic::swapWavefields)
      .def(
          "get_current_field",
          [](const SEMsolverDataElastic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const SEMsolverDataElastic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPreviousField(i));
          },
          py::arg("i"))
      .def(
          "get_prevprev_field",
          [](const SEMsolverDataElastic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPrevPrevField(i));
          },
          py::arg("i"))
      .def("print", &SEMsolverDataElastic::print);
}

void bind_acoustoelastic_solver_data(py::module_ &m) {
  py::class_<SEMsolverDataAcoustoElastic, Solver::DataStruct, std::shared_ptr<SEMsolverDataAcoustoElastic>>(
      m, "SEMsolverDataAcoustoElastic")
      .def(py::init<const WavefieldAcoustoElastic &, const RhsAcoustoElastic &>(), py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &SEMsolverDataAcoustoElastic::swapWavefields)
      .def(
          "get_current_field",
          [](const SEMsolverDataAcoustoElastic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const SEMsolverDataAcoustoElastic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPreviousField(i));
          },
          py::arg("i"))
      .def(
          "get_prevprev_field",
          [](const SEMsolverDataAcoustoElastic &self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPrevPrevField(i));
          },
          py::arg("i"))
      .def("print", &SEMsolverDataAcoustoElastic::print);
}

}  // namespace solver
#endif  // FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
