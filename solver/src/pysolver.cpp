#include <pybind11/pybind11.h>
#include "SEMsolver.hpp"

namespace py = pybind11;
PYBIND11_MODULE(pysolver, m) {
  py::class_<SEMsolver>(m, "SEMSolver")
    .def("computeOneStep", &SEMsolver::computeOneStep)
    .def("outputPnValues", &SEMsolver::outputPnValues);
}
