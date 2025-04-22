#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SEMproxy.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pysem, m) {
    m.def("initialize_kokkos", []() {
        if (!Kokkos::is_initialized()) {
            Kokkos::initialize();
        }
    });

    m.def("finalize_kokkos", []() {
        if (Kokkos::is_initialized()) {
            Kokkos::finalize();
        }
    });

    py::class_<SEMproxy>(m, "SEMproxy")
        .def(py::init<int, int, int, float>())
        .def("get_mesh_info", &SEMproxy::getMeshInfo)
        .def("run", &SEMproxy::run)
        .def("initFiniteElem", &SEMproxy::initFiniteElem);
}
