#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <commonMacros.hpp>
#include <SEMsolver.hpp>
#include <SolverBase.hpp>
#include <solverFactory.hpp>

namespace py = pybind11;

PYBIND11_MODULE(solver, m) {

  // Create submodule 'solver'
  m.attr("__name__") = "pyproxys.solver";

  // Bind enums
  py::enum_<SolverFactory::methodType>(m, "MethodType")
    .value("SEM", SolverFactory::SEM)
    .value("DG", SolverFactory::DG);

  py::enum_<SolverFactory::implemType>(m, "ImplemType")
    .value("CLASSIC", SolverFactory::CLASSIC)
    .value("GEOS", SolverFactory::GEOS)
    .value("OPTIM", SolverFactory::OPTIM)
    .value("SHIVA", SolverFactory::SHIVA);

  py::enum_<SolverFactory::meshType>(m, "MeshType")
    .value("STRUCT", SolverFactory::Struct)
    .value("UNSTRUCT", SolverFactory::Unstruct);

  // Bind DataStruct
  py::class_<SolverBase::DataStruct, std::shared_ptr<SolverBase::DataStruct>>(m, "DataStruct")
    .def(py::init<>());

  // Bind SEMsolverData (inherits from DataStruct)
  py::class_<SEMsolverData, SolverBase::DataStruct, std::shared_ptr<SEMsolverData>>(m, "SEMsolverData")
    .def(py::init<int, int,
                  const ARRAY_REAL_VIEW&,
                  const ARRAY_REAL_VIEW&,
                  const VECTOR_INT_VIEW&,
                  ARRAY_REAL_VIEW&>(),
         py::arg("i1"),
         py::arg("i2"),
         py::arg("rhs_term"),
         py::arg("pn_global"),
         py::arg("rhs_element"),
         py::arg("rhs_weights"))
    .def_readwrite("i1", &SEMsolverData::m_i1)
    .def_readwrite("i2", &SEMsolverData::m_i2)
    .def_property_readonly("rhs_term", [](SEMsolverData &self) { return self.m_rhsTerm; })
    .def_property_readonly("pn_global", [](SEMsolverData &self) { return self.m_pnGlobal; })
    .def_property_readonly("rhs_element", [](SEMsolverData &self) { return self.m_rhsElement; })
    .def_property_readonly("rhs_weights", [](SEMsolverData &self) { return self.m_rhsWeights; });

  // Bind SolverBase
  py::class_<SolverBase, std::shared_ptr<SolverBase>>(m, "SolverBase")
    .def("compute_fe_init", &SolverBase::computeFEInit,
        py::arg("model"))
    .def("compute_one_step", &SolverBase::computeOneStep,
        py::arg("dt"),
        py::arg("time_sample"),
        py::arg("data"))
    .def("output_pn_values", &SolverBase::outputPnValues,
        py::arg("index_time_step"),
        py::arg("i1"),
        py::arg("my_element_source"),
        py::arg("pn_global"));

  // Bind Solver factory function (returns shared_ptr<SolverBase>)
  m.def("create_solver",
    [](SolverFactory::methodType methodType,
       SolverFactory::implemType implemType,
       SolverFactory::meshType meshType,
       int order) {
        auto solver = SolverFactory::createSolver(methodType, implemType, meshType, order);
        return std::shared_ptr<SolverBase>(std::move(solver)); // pyfwi needs to do solver2 = solver1
    },
    py::arg("method_type"),
    py::arg("implem_type"),
    py::arg("mesh_type"),
    py::arg("order"));

}
