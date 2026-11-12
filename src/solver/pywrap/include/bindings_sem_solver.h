#ifndef FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <memory>

#include "common_macros.h"
#include "field_sampler.h"
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

void bind_field_sampler(py::module_ &m) {
  py::class_<FieldSampler, std::shared_ptr<FieldSampler>>(m, "FieldSampler",
                                                         "Weighted sums of a nodal field at fixed points, on the device "
                                                         "(receivers): only the point values are copied to the host.")
      .def(py::init([](const py::array_t<int, py::array::c_style | py::array::forcecast> &nodes,
                       const py::array_t<float, py::array::c_style | py::array::forcecast> &weights,
                       int nodes_per_point) {
             if (nodes.size() != weights.size())
               throw std::invalid_argument("FieldSampler: nodes and weights differ in size");
             vectorInt vn = allocateVector<vectorInt>(nodes.size(), "field_sampler_nodes");
             vectorReal vw = allocateVector<vectorReal>(weights.size(), "field_sampler_weights");
             Kokkos::deep_copy(vn, Kokkos::View<const int *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                                       nodes.data(), nodes.size()));
             Kokkos::deep_copy(vw, Kokkos::View<const float *, Kokkos::HostSpace,
                                                Kokkos::MemoryTraits<Kokkos::Unmanaged>>(weights.data(), weights.size()));
             return std::make_shared<FieldSampler>(vn, vw, nodes_per_point);
           }),
           py::arg("nodes"), py::arg("weights"), py::arg("nodes_per_point"))
      .def_property_readonly("num_points", &FieldSampler::numPoints)
      .def(
          "sample",
          [](const FieldSampler &self, const SEMsolverDataAcoustic &data, int component) {
            self.sample(data.getCurrentField(component));
            py::array_t<float> out(self.numPoints());
            Kokkos::deep_copy(Kokkos::View<float *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                                  out.mutable_data(), self.numPoints()),
                              self.values());
            return out;
          },
          py::arg("data"), py::arg("component") = 0, "Values at the points of the current field of an acoustic SEM data.");
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
