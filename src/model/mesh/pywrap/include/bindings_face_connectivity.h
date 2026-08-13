#ifndef FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_FACE_CONNECTIVITY_H_
#define FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_FACE_CONNECTIVITY_H_
#include <pybind11/pybind11.h>

#include <KokkosExp_InterOp.hpp>

#include "face_connectivity_unstruct.h"
#include "model_unstruct.h"

namespace py = pybind11;

namespace bindings {

template <typename FloatType, typename ScalarType>
void bindFaceConnectivityUnstruct(py::module &m) {
  using FaceConnData = model::FaceConnectivityUnstructData<FloatType, ScalarType>;
  using FaceConn = model::FaceConnectivityUnstruct<FloatType, ScalarType>;
  using MeshType = model::ModelUnstruct<FloatType, ScalarType>;

  std::string data_name = model::model_class_name<FloatType, ScalarType>("FaceConnectivityUnstructData");
  std::string class_name = model::model_class_name<FloatType, ScalarType>("FaceConnectivityUnstruct");

  // Bind Data struct
  py::class_<FaceConnData>(m, data_name.c_str())
      .def(py::init<>())
      .def_readwrite("n_faces", &FaceConnData::n_faces)
      .def_readwrite("ndofs_per_face", &FaceConnData::ndofs_per_face)
      .def_property(
          "elem_to_faces",
          [](FaceConnData &self) -> Kokkos::Experimental::python_view_type_t<decltype(self.elem_to_faces)> {
            return self.elem_to_faces;
          },
          [](FaceConnData &self, Kokkos::Experimental::python_view_type_t<decltype(self.elem_to_faces)> v) {
            self.elem_to_faces = v;
          })

      .def_property(
          "face_dofs",
          [](FaceConnData &self) -> Kokkos::Experimental::python_view_type_t<decltype(self.face_dofs)> {
            return self.face_dofs;
          },
          [](FaceConnData &self, Kokkos::Experimental::python_view_type_t<decltype(self.face_dofs)> v) {
            self.face_dofs = v;
          })

      .def_property(
          "face_elem_owner",
          [](FaceConnData &self) -> Kokkos::Experimental::python_view_type_t<decltype(self.face_elem_owner)> {
            return self.face_elem_owner;
          },
          [](FaceConnData &self, Kokkos::Experimental::python_view_type_t<decltype(self.face_elem_owner)> v) {
            self.face_elem_owner = v;
          })

      .def_property(
          "face_elem_neighbor",
          [](FaceConnData &self) -> Kokkos::Experimental::python_view_type_t<decltype(self.face_elem_neighbor)> {
            return self.face_elem_neighbor;
          },
          [](FaceConnData &self, Kokkos::Experimental::python_view_type_t<decltype(self.face_elem_neighbor)> v) {
            self.face_elem_neighbor = v;
          })

      .def_property(
          "face_local_owner",
          [](FaceConnData &self) -> Kokkos::Experimental::python_view_type_t<decltype(self.face_local_owner)> {
            return self.face_local_owner;
          },
          [](FaceConnData &self, Kokkos::Experimental::python_view_type_t<decltype(self.face_local_owner)> v) {
            self.face_local_owner = v;
          })

      .def_property(
          "face_local_neighbor",
          [](FaceConnData &self) -> Kokkos::Experimental::python_view_type_t<decltype(self.face_local_neighbor)> {
            return self.face_local_neighbor;
          },
          [](FaceConnData &self, Kokkos::Experimental::python_view_type_t<decltype(self.face_local_neighbor)> v) {
            self.face_local_neighbor = v;
          });
  // Bind Class
  py::class_<FaceConn>(m, class_name.c_str())
      .def(py::init<>())
      .def(py::init<const FaceConnData &>())
      .def("build", &FaceConn::template build<MeshType>)
      .def("get_number_of_faces", &FaceConn::getNumberOfFaces)
      .def("is_boundary_face", &FaceConn::isBoundaryFace)
      .def("get_global_face", &FaceConn::getGlobalFace)
      .def("get_global_node_from_face", &FaceConn::getGlobalNodeFromFace);
}

}  // namespace bindings
#endif  // FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_FACE_CONNECTIVITY_H_