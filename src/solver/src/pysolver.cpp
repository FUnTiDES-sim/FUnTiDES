#include "SEMsolver.hpp"
#include <libpykokkos.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// void compute_debug_args(py::object self_obj, py::object a0, py::object a1,
//                         py::object a2, py::object a3, py::object a4,
//                         py::object a5, py::object a6, py::object a7,
//                         py::object a8) {
//   auto &solver = self_obj.cast<SEMsolver &>();

//   int t = a0.cast<int>();
//   int order = a1.cast<int>();
//   int npts = a2.cast<int>();
//   int i1 = a3.cast<int>();
//   int i2 = a4.cast<int>();
//   SEMinfo &info = a5.cast<SEMinfo &>();
//   auto kk1 = a6.cast<
//       Kokkos::View<float **, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace>>();
//   auto kk2 = a7.cast<
//       Kokkos::View<float **, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace>>();
//   auto kk3 =
//       a8.cast<Kokkos::View<int *, Kokkos::LayoutLeft,
//       Kokkos::CudaUVMSpace>>();

//   solver.computeOneStep(t, order, npts, i1, i2, info, kk1, kk2, kk3);
// }

PYBIND11_MODULE(pysolver, m) {

  py::object pysem = py::module_::import("pysem");
  py::object SEMinfoClass = pysem.attr("SEMinfo");

  py::class_<SEMsolver>(m, "SEMsolver")
      .def(py::init<>())

      // Basic getters/setters
      .def("get_order", &SEMsolver::getOrder)
      .def("set_order", &SEMsolver::setOrder)
      .def("get_vmin", &SEMsolver::getVMin)
      .def("set_vmin", &SEMsolver::setVMin)

      // Kokkos::View container accessors
      .def("get_model", &SEMsolver::getModel)
      .def("set_model", &SEMsolver::setModel_wrapper)
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

      // For internal mesh
      .def("computeFEInit",&SEMsolver::computeFEInit)

      // For pyFWI
      .def("computeFEInitWithoutMesh",    &SEMsolver::computeFEInitWithoutMesh_) // = allocate + init
      .def("allocateFEarraysWithoutMesh", &SEMsolver::allocateFEarraysWithoutMesh)
      .def("initFEarraysWithoutMesh",     &SEMsolver::initFEarraysWithoutMesh_)

      .def("computeOneStep", &SEMsolver::computeOneStep_wrapper);

  // m.def("compute_debug_args", &compute_debug_args);

  // .def("compute_one_step", [](SEMsolver &self, int t, int order,
  //                             int nPointsPerElement, int i1, int i2,
  //                             SEMinfo &info, py::array_t<float> rhsTerm,
  //                             py::array_t<float> pnGlobal,
  //                             py::array_t<int> rhsElement) {
  //   // Convert numpy arrays to raw pointers and wrap in Kokkos::View
  //   auto buf1 = rhsTerm.request();
  //   auto buf2 = pnGlobal.request();
  //   auto buf3 = rhsElement.request();

  //   // Créer et remplir des vues dans CudaUVMSpace
  //   Kokkos::View<float **, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace> rhs(
  //       "rhs", buf1.shape[0], buf1.shape[1]);
  //   Kokkos::View<float **, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace> pn(
  //       "pn", buf2.shape[0], buf2.shape[1]);
  //   Kokkos::View<int *, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace> elem(
  //       "elem", buf3.shape[0]);

  //   // Créer vues temporaires HostSpace et copier les données depuis
  //   // NumPy
  //   Kokkos::View<float **, Kokkos::LayoutLeft, Kokkos::HostSpace> rhs_host(
  //       static_cast<float *>(buf1.ptr), buf1.shape[0], buf1.shape[1]);
  //   Kokkos::View<float **, Kokkos::LayoutLeft, Kokkos::HostSpace> pn_host(
  //       static_cast<float *>(buf2.ptr), buf2.shape[0], buf2.shape[1]);
  //   Kokkos::View<int *, Kokkos::LayoutLeft, Kokkos::HostSpace> elem_host(
  //       static_cast<int *>(buf3.ptr), buf3.shape[0]);

  //   Kokkos::deep_copy(rhs, rhs_host);
  //   Kokkos::deep_copy(pn, pn_host);
  //   Kokkos::deep_copy(elem, elem_host);

  //   // Appel de la méthode native
  //   Kokkos::fence("before computeOneStep");
  //   self.computeOneStep(t, order, nPointsPerElement, i1, i2, info, rhs, pn,
  //                       elem);
  //   Kokkos::fence("after computeOneStep");

  //   Kokkos::deep_copy(rhs_host, rhs);
  //   Kokkos::deep_copy(pn_host, pn);
  //   Kokkos::deep_copy(elem_host, elem);

  //   Kokkos::fence();
  // });
}
