#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <memory>

#include "common_macros.h"
#include "rhs_acoustic.h"
#include "rhs_elastic.h"
#include "sem_solver.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace py = pybind11;

namespace solver
{
namespace fe
{

void bind_data_struct(py::module_ &m)
{
  py::class_<SolverBase::DataStruct, std::shared_ptr<SolverBase::DataStruct>>(
      m, "DataStruct")
      .def("print", &SolverBase::DataStruct::print);
}

void bind_acoustic_solver_data(py::module_ &m)
{
  py::class_<SEMsolverDataAcoustic, SolverBase::DataStruct,
             std::shared_ptr<SEMsolverDataAcoustic>>(m, "SEMsolverDataAcoustic")
      .def(py::init([](std::shared_ptr<WavefieldAcoustic> wavefield,
                       std::shared_ptr<RhsAcoustic> rhs) {
             // Get raw pointers from shared_ptr - Python manages lifetime
             return new SEMsolverDataAcoustic(wavefield.get(), rhs.get());
           }),
           py::arg("wavefield"), py::arg("rhs"),
           py::keep_alive<1, 2>(),  // Keep wavefield alive with SEMsolverDataAcoustic
           py::keep_alive<1, 3>())  // Keep rhs alive with SEMsolverDataAcoustic
      .def("print", &SEMsolverDataAcoustic::print);
}

void bind_elastic_solver_data(py::module_ &m)
{
  py::class_<SEMsolverDataElastic, SolverBase::DataStruct,
             std::shared_ptr<SEMsolverDataElastic>>(m, "SEMsolverDataElastic")
      .def(py::init([](std::shared_ptr<WavefieldElastic> wavefield,
                       std::shared_ptr<RhsElastic> rhs) {
             // Get raw pointers from shared_ptr - Python manages lifetime
             return new SEMsolverDataElastic(wavefield.get(), rhs.get());
           }),
           py::arg("wavefield"), py::arg("rhs"),
           py::keep_alive<1, 2>(),  // Keep wavefield alive with SEMsolverDataElastic
           py::keep_alive<1, 3>())  // Keep rhs alive with SEMsolverDataElastic
      .def("print", &SEMsolverDataElastic::print);
}

void bind_sem_solver_base(py::module_ &m)
{
  py::class_<SEMSolverBase, std::shared_ptr<SEMSolverBase>>(m, "SEMSolverBase")
      .def("compute_fe_init", &SEMSolverBase::computeFEInit, py::arg("model"),
           py::arg("sponge_size") = std::array<float, 3>{0.0f, 0.0f, 0.0f},
           py::arg("sponge_surface") = true, py::arg("taper_delta") = 0)
      .def("compute_one_step", &SEMSolverBase::computeOneStep, py::arg("dt"),
           py::arg("time_sample"), py::arg("data"))
      .def("output_solution_values", &SEMSolverBase::outputSolutionValues,
           py::arg("t"), py::arg("e"), py::arg("field_global"),
           py::arg("field_name"));
}

void bind_solver_factory(py::module_ &m)
{
  m.def(
      "create_solver",
      [](methodType method, implemType implem, meshType mesh,
         modelLocationType modelLocation, physicType physic, int order) {
        auto solver = SolverFactory::createSolver(method, implem, mesh,
                                                  modelLocation, physic, order);
        return std::shared_ptr<SEMSolverBase>(std::move(solver));
      },
      py::arg("method_type"), py::arg("implem_type"), py::arg("mesh_type"),
      py::arg("model_location"), py::arg("physic_type"), py::arg("order"));
}

}  // namespace fe
}  // namespace solver
