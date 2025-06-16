#include "SEMsolver.hpp"
#include <libpykokkos.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pysolver, m) {

  py::class_<SEMsolver>(m, "SEMsolver")
      .def(py::init<>())

      // Basic getters/setters
      .def("get_order", &SEMsolver::getOrder)
      .def("set_order", &SEMsolver::setOrder)
      .def("get_vmin", &SEMsolver::getVMin)
      .def("set_vmin", &SEMsolver::setVMin)

      // Kokkos::View container accessors
      .def("get_model", &SEMsolver::getModel)
      .def("set_model", &SEMsolver::setModel)
      .def("get_mass_matrix_global", &SEMsolver::getMassMatrixGlobal)
      .def("set_mass_matrix_global", &SEMsolver::setMassMatrixGlobal)
      .def("get_global_nodes_list", &SEMsolver::getGlobalNodesList)
      .def("set_global_nodes_list", &SEMsolver::setGlobalNodesList)
      .def("get_list_of_interior_nodes", &SEMsolver::getListOfInteriorNodes)
      .def("set_list_of_interior_nodes", &SEMsolver::setListOfInteriorNodes)
      .def("get_list_of_damping_nodes", &SEMsolver::getListOfDampingNodes)
      .def("set_list_of_damping_nodes", &SEMsolver::setListOfDampingNodes)
      .def("get_sponge_taper_coeff", &SEMsolver::getSpongeTaperCoeff)
      .def("set_sponge_taper_coeff", &SEMsolver::setSpongeTaperCoeff)

      .def("initFEarrays", &SEMsolver::initFEarrays)
      .def("compute_one_step", &SEMsolver::computeOneStep);
}
