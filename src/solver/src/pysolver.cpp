#include "KokkosExp_InterOp.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "SEMsolver.hpp"
#include "dataType.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using SEMmesh = CartesianSEMmesh<float, float, int, int, 2>;

PYBIND11_MODULE(pysolver, m) {

  py::class_<SEMmesh>(m, "SEMmesh")
      .def(py::init<>())
      .def("get_nb_points_per_element", &SEMmesh::getNumberOfPointsPerElement)
      .def(py::init<int, int, int, float, float, float, int>());

  py::class_<SEMsolver>(m, "SEMsolver")
      .def(py::init<>())
      .def(py::init<const SEMmesh &>())

      .def("computeFEInit", &SEMsolver::computeFEInit)

      .def(
          "computeOneStep",
          [](SEMsolver &self, int t, int i1, int i2,
             Kokkos::Experimental::python_view_type_t<
                 Kokkos::View<float **, Layout, MemSpace>>
                 rhsTerm_view,
             Kokkos::Experimental::python_view_type_t<
                 Kokkos::View<float **, Layout, MemSpace>>
                 pnGlobal_view,
             Kokkos::Experimental::python_view_type_t<
                 Kokkos::View<int *, Layout, MemSpace>>
                 rhsElement_view,
             Kokkos::Experimental::python_view_type_t<Kokkos::View<float **, Layout, MemSpace>> rhsWeights) {
            self.computeOneStep(t, i1, i2, rhsTerm_view, pnGlobal_view,
                                rhsElement_view, rhsWeights);
          },
          "Compute one step of the solver");
}
