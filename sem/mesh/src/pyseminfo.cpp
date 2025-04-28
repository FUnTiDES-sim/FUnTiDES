#include <pybind11/pybind11.hpp>
#include "SEMdata.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyseminfo, m) {
  py::class_<SEMinfo>(m, "SEMinfo")
        .def(py::init<>())
        .def_readwrite("numberOfNodes", &SEMinfo::numberOfNodes)
        .def_readwrite("numberOfElements", &SEMinfo::numberOfElements)
        .def_readwrite("numberOfPointsPerElement", &SEMinfo::numberOfPointsPerElement)
        .def_readwrite("numberOfInteriorNodes", &SEMinfo::numberOfInteriorNodes)
        .def_readonly("myNumberOfRHS", &SEMinfo::myNumberOfRHS)
        .def_readonly_static("myOrderNumber", []() { return SEMinfo::myOrderNumber; })
        .def_readonly("myTimeStep", &SEMinfo::myTimeStep)
        .def_property_readonly("nPointsPerElement", [](const SEMinfo &self) { return self.nPointsPerElement; })
        .def_readonly("f0", &SEMinfo::f0)
        .def_readonly("myTimeMax", &SEMinfo::myTimeMax)
        .def_readonly("sourceOrder", &SEMinfo::sourceOrder)
        .def_property_readonly("myNumSamples", [](const SEMinfo &self) { return self.myNumSamples; })
        .def_readwrite("myElementSource", &SEMinfo::myElementSource);
}
