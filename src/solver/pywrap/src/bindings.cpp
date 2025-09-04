#include <KokkosExp_InterOp.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dataType.hpp>
#include <SEMsolver.hpp>
#include <SolverBase.hpp>

#include "solverFactory.hpp"

namespace py = pybind11;

// Define type aliases for Kokkos views
using Float2DView = Kokkos::View<float **, Layout, MemSpace>;
using Int1DView = Kokkos::View<int *, Layout, MemSpace>;
using PyFloat2DView = Kokkos::Experimental::python_view_type_t<Float2DView>;
using PyInt1DView = Kokkos::Experimental::python_view_type_t<Int1DView>;

PYBIND11_MODULE(pyproxys, m) {

  // Bind enums
  py::enum_<SolverFactory::methodType>(m, "MethodType")
    .value("SEM", SolverFactory::SEM)
    .value("DG", SolverFactory::DG)
    .export_values();

  py::enum_<SolverFactory::implemType>(m, "ImplemType")
    .value("CLASSIC", SolverFactory::CLASSIC)
    .value("GEOS", SolverFactory::GEOS)
    .value("OPTIM", SolverFactory::OPTIM)
    .value("SHIVA", SolverFactory::SHIVA)
    .export_values();

  py::enum_<SolverFactory::meshType>(m, "MeshType")
    .value("Struct", SolverFactory::Struct)
    .value("Unstruct", SolverFactory::Unstruct)
    .export_values();

  // Bind SolverBase

    py::class_<SolverBase::DataStruct, std::shared_ptr<SolverBase::DataStruct>>(m, "DataStruct")
      .def(py::init<>());

    py::class_<SolverBase, std::shared_ptr<SolverBase>>(m, "SolverBase")
      .def("computeFEInit", &SolverBase::computeFEInit)
      .def("computeOneStep",
        [](SolverBase &self, float dt, int timeSample, SolverBase::DataStruct &data) {
          self.computeOneStep(dt, timeSample, data);
        }
      )
      .def("outputPnValues", &SolverBase::outputPnValues);

  // Bind factory function
  m.def("createSolver",
    [](SolverFactory::methodType methodType,
       SolverFactory::implemType implemType,
       SolverFactory::meshType meshType,
       int order) {
        auto solver = SolverFactory::createSolver(methodType, implemType, meshType, order);
        return std::shared_ptr<SolverBase>(std::move(solver));
    },
    py::arg("methodType"),
    py::arg("implemType"),
    py::arg("meshType"),
    py::arg("order"),
    "Create a solver instance using the factory."
  );

}
