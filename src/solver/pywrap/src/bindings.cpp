#include <common_macros.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sem_solver_acoustic.h>
#include <sem_solver_elastic.h>
#include <solver_factory.h>

#include <KokkosExp_InterOp.hpp>

namespace py = pybind11;

PYBIND11_MODULE(solver, m)
{
  // Create submodule 'solver'
  m.attr("__name__") = "pyfuntides.solver";

  // Bind enums
  // (!they are not python enums, just to select the template at runtime!)
  py::enum_<SolverFactory::methodType>(m, "MethodType")
      .value("SEM", SolverFactory::SEM)
      .value("DG", SolverFactory::DG);

  py::enum_<SolverFactory::implemType>(m, "ImplemType")
      .value("MAKUTU", SolverFactory::MAKUTU)
      .value("SHIVA", SolverFactory::SHIVA);

  py::enum_<SolverFactory::meshType>(m, "MeshType")
      .value("STRUCT", SolverFactory::Struct)
      .value("UNSTRUCT", SolverFactory::Unstruct);

  py::enum_<SolverFactory::modelLocationType>(m, "ModelLocationType")
      .value("ONNODES", SolverFactory::modelLocationType::OnNodes)
      .value("ONELEMENTS", SolverFactory::modelLocationType::OnElements)
      .export_values();

  py::enum_<SolverFactory::physicType>(m, "PhysicType")
      .value("ACOUSTIC", SolverFactory::Acoustic)
      .value("ELASTIC", SolverFactory::Elastic)
      .export_values();

  // Bind DataStruct
  py::class_<SolverBase::DataStruct, std::shared_ptr<SolverBase::DataStruct>>(
      m, "DataStruct")
      .def("print", &SolverBase::DataStruct::print);

  // Bind SEMsolverDataAcoustic (inherits from SolverBase::DataStruct)
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

  // Bind SEMSolverDataElastic (inherits from SolverBase::DataStruct)
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

  // Bind SEMSolverBase
  py::class_<SEMSolverBase, std::shared_ptr<SEMSolverBase>>(m, "SEMSolverBase")
      .def("compute_fe_init", &SEMSolverBase::computeFEInit, py::arg("model"),
           py::arg("sponge_size") = std::array<float, 3>{0.0f, 0.0f, 0.0f},
           py::arg("sponge_surface") = true, py::arg("taper_delta") = 0)
      .def("compute_one_step", &SEMSolverBase::computeOneStep, py::arg("dt"),
           py::arg("time_sample"), py::arg("data"))
      .def("output_solution_values", &SEMSolverBase::outputSolutionValues,
           py::arg("index_time_step"), py::arg("i1"),
           py::arg("my_element_source"), py::arg("field_global"),
           py::arg("field_name"));

  // Bind Solver factory function (returns shared_ptr<SolverBase>)
  m.def(
      "create_solver",
      [](SolverFactory::methodType methodType,
         SolverFactory::implemType implemType, SolverFactory::meshType meshType,
         SolverFactory::modelLocationType modelLocation,
         SolverFactory::physicType physicType, int order) {
        auto solver = SolverFactory::createSolver(
            methodType, implemType, meshType, modelLocation, physicType, order);
        return std::shared_ptr<SEMSolverBase>(
            std::move(solver));  // pyfwi needs to do solver2 = solver1
      },
      py::arg("method_type"), py::arg("implem_type"), py::arg("mesh_type"),
      py::arg("model_location"), py::arg("physic_type"), py::arg("order"));
}
