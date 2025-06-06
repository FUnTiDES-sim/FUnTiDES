#include "SEMsolver.hpp"
#include <libpykokkos.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pk = pykokkos;

PYBIND11_MODULE(sem_bindings, m) {
  pk::initialize(); // Initializes Kokkos

  py::class_<SEMsolver>(m, "SEMsolver")
      .def(py::init<>())
      .def("computeFEInit", [](SEMsolver &self, SEMinfo &info,
                               Mesh &mesh) { self.computeFEInit(info, mesh); })
      .def("computeOneStep",
           [](SEMsolver &self, int timeSample, int order, int nPointsPerElement,
              int i1, int i2, SEMinfo &info, pk::view_type<float **> rhsTerm,
              pk::view_type<float **> pnGlobal, std::vector<int> rhsElement) {
             self.computeOneStep(timeSample, order, nPointsPerElement, i1, i2,
                                 info, rhsTerm, pnGlobal, rhsElement);
           });

  // Register pykokkos view types (arrayReal2D, etc.)
  pk::add_view_to_module<float, 2>(m, "arrayReal2D");
  pk::add_view_to_module<float, 2>(m, "pnGlobal");
  pk::add_view_to_module<float, 2>(m, "rhsTerm");
}
