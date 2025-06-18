#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SEMproxy.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pysem, m) {
  py::class_<SEMproxy>(m, "SEMproxy")
      .def(py::init<int, int, int, float>())
      .def("get_mesh_info", &SEMproxy::getMeshInfo)
      .def("run", &SEMproxy::run)
      .def("initFiniteElem", &SEMproxy::initFiniteElem)
      .def("get_my_rhs_term", &SEMproxy::getMyRHSTerm)
      .def("set_my_rhs_term", &SEMproxy::setMyRHSTerm)
      .def("get_pn_global", &SEMproxy::getPnGlobal)
      .def("set_pn_global", &SEMproxy::setPnGlobal)
      .def("get_rhs_element", &SEMproxy::getRhsElement)
      .def("set_rhs_element", &SEMproxy::setRhsElement);

  py::class_<SEMmesh>(m, "SEMmesh")
      .def(py::init<>())
      .def(py::init<int, int, int, float, float, float, int, int, bool>())
      .def("getSpongeSize", &SEMmesh::getSpongeSize)
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

  py::class_<SEMinfo>(m, "SEMinfo")
      .def(py::init<>()) // you'll need a default constructor or custom init
      .def_readwrite("numberOfNodes", &SEMinfo::numberOfNodes)
      .def_readwrite("numberOfElements", &SEMinfo::numberOfElements)
      .def_readwrite("numberOfPointsPerElement",
                     &SEMinfo::numberOfPointsPerElement)
      .def_readwrite("numberOfInteriorNodes", &SEMinfo::numberOfInteriorNodes)
      .def_readwrite("numberOfDampingNodes", &SEMinfo::numberOfDampingNodes)
      .def_readwrite("numberOfSpongeNodes", &SEMinfo::numberOfSpongeNodes)
      .def_readwrite("myNumSamples", &SEMinfo::myNumSamples)
      .def_readwrite("myElementSource", &SEMinfo::myElementSource);
}
