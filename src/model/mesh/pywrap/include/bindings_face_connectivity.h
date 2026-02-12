#ifndef BINDINGS_FACE_CONNECTIVITY_H_
#define BINDINGS_FACE_CONNECTIVITY_H_

#include <pybind11/pybind11.h>

#include "face_connectivity_unstruct.h"

namespace py = pybind11;

namespace bindings
{

template <typename FloatType, typename ScalarType>
void bindFaceConnectivityUnstruct(py::module& m)
{
  using FaceConnData =
      model::FaceConnectivityUnstructData<FloatType, ScalarType>;
  using FaceConn = model::FaceConnectivityUnstruct<FloatType, ScalarType>;

  std::string data_name = model::model_class_name<FloatType, ScalarType>(
      "FaceConnectivityUnstructData");
  std::string class_name = model::model_class_name<FloatType, ScalarType>(
      "FaceConnectivityUnstruct");

  // Bind Data struct
  py::class_<FaceConnData>(m, data_name.c_str())
      .def(py::init<>())
      .def_readwrite("n_faces", &FaceConnData::n_faces)
      .def_readwrite("ndofs_per_face", &FaceConnData::ndofs_per_face)
      .def_readwrite("elem_to_faces", &FaceConnData::elem_to_faces)
      .def_readwrite("face_dofs", &FaceConnData::face_dofs)
      .def_readwrite("face_elem_owner", &FaceConnData::face_elem_owner)
      .def_readwrite("face_elem_neighbor", &FaceConnData::face_elem_neighbor)
      .def_readwrite("face_local_owner", &FaceConnData::face_local_owner)
      .def_readwrite("face_local_neighbor", &FaceConnData::face_local_neighbor);

  // Bind Class
  py::class_<FaceConn>(m, class_name.c_str())
      .def(py::init<>())
      .def(py::init<const FaceConnData&>())
      .def("build", &FaceConn::build)
      .def("get_number_of_faces", &FaceConn::getNumberOfFaces)
      .def("is_boundary_face", &FaceConn::isBoundaryFace)
      .def("get_global_face", &FaceConn::getGlobalFace)
      .def("get_global_node_from_face", &FaceConn::getGlobalNodeFromFace);
}

}  // namespace bindings
#endif  // BINDINGS_FACE_CONNECTIVITY_H_