#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "common_macros.h"
#include "sem_solver.h"
#include "solver_factory.h"

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
      .def(
          py::init<int, int,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<VECTOR_INT_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>>(),
          py::arg("i1"), py::arg("i2"), py::arg("rhs_term"),
          py::arg("pn_global"), py::arg("rhs_element"), py::arg("rhs_weights"))
      .def("print", &SEMsolverDataAcoustic::print)
      .def_readwrite("i1", &SEMsolverDataAcoustic::m_i1)
      .def_readwrite("i2", &SEMsolverDataAcoustic::m_i2);
}

void bind_elastic_solver_data(py::module_ &m)
{
  py::class_<SEMsolverDataElastic, SolverBase::DataStruct,
             std::shared_ptr<SEMsolverDataElastic>>(m, "SEMsolverDataElastic")
      .def(
          py::init<int, int,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>,
                   Kokkos::Experimental::python_view_type_t<VECTOR_INT_VIEW>,
                   Kokkos::Experimental::python_view_type_t<ARRAY_REAL_VIEW>>(),
          py::arg("i1"), py::arg("i2"), py::arg("rhs_termx"),
          py::arg("rhs_termy"), py::arg("rhs_termz"), py::arg("uxn_global"),
          py::arg("uyn_global"), py::arg("uzn_global"), py::arg("rhs_element"),
          py::arg("rhs_weights"))
      .def("print", &SEMsolverDataElastic::print)
      .def_readwrite("i1", &SEMsolverDataElastic::m_i1)
      .def_readwrite("i2", &SEMsolverDataElastic::m_i2);
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
