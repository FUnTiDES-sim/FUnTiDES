#include "KokkosExp_InterOp.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "SEMsolver.hpp"
#include "dataType.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using SEMmesh = cartesianSEMmesh<float, int, int, 2>;

PYBIND11_MODULE(pysolver, m) {

  py::class_<SEMmesh>(m, "SEMmesh")
      .def(py::init<>())
      .def(py::init<int, int, int, float, float, float, int, int, bool>())
      .def("getNumberOfNodes", &SEMmesh::getNumberOfNodes)
      .def("getNumberOfElements", &SEMmesh::getNumberOfElements)
      .def("getNumberOfInteriorNodes",
           py::overload_cast<>(&SEMmesh::getNumberOfInteriorNodes, py::const_))
      .def("getNumberOfInteriorNodes",
           py::overload_cast<int>(&SEMmesh::getNumberOfInteriorNodes,
                                  py::const_),
           py::arg("spongeSize"))
      .def("getNumberOfPointsPerElement",
           &SEMmesh::getNumberOfPointsPerElement);

  // py::class_<SEMinfo>(m, "SEMinfo")
  //     .def(py::init<>()) // you'll need a default constructor or custom init
  //     .def_readwrite("numberOfNodes", &SEMinfo::numberOfNodes)
  //     .def_readwrite("numberOfElements", &SEMinfo::numberOfElements)
  //     .def_readwrite("numberOfPointsPerElement",
  //                    &SEMinfo::numberOfPointsPerElement)
  //     .def_readwrite("numberOfInteriorNodes",
  //     &SEMinfo::numberOfInteriorNodes) .def_readwrite("numberOfSpongeNodes",
  //     &SEMinfo::numberOfSpongeNodes) .def_readwrite("myNumSamples",
  //     &SEMinfo::myNumSamples) .def_readwrite("myElementSource",
  //     &SEMinfo::myElementSource);

  py::class_<SEMsolver>(m, "SEMsolver")
      .def(py::init<>())

      // Kokkos::View container accessors
      .def("set_model",
           [](SEMsolver &self,
              py::array_t<float, py::array::c_style | py::array::forcecast>
                  model_np) {
             auto model_view =
                 Kokkos::View<float *, Kokkos::LayoutRight, MemSpace>(
                     static_cast<float *>(model_np.mutable_data()),
                     model_np.shape(0));
             self.setModel(model_view);
             Kokkos::fence();
           })

      .def("computeFEInit", &SEMsolver::computeFEInit)

      .def(
          "computeOneStep",
          [](SEMsolver &self, int t, int order, int npts, int i1, int i2,
             SEMinfo info,
             Kokkos::Experimental::python_view_type_t<
                 Kokkos::View<float **, Layout, MemSpace>>
                 rhsTerm_view,
             Kokkos::Experimental::python_view_type_t<
                 Kokkos::View<float **, Layout, MemSpace>>
                 pnGlobal_view,
             Kokkos::Experimental::python_view_type_t<
                 Kokkos::View<int *, Layout, MemSpace>>
                 rhsElement_view) {
            self.computeOneStep(t, order, npts, i1, i2, info, rhsTerm_view,
                                pnGlobal_view, rhsElement_view);
          },
          "Compute one step of the solver");
}
