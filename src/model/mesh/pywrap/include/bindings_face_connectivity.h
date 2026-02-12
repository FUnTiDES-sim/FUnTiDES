#ifndef BINDINGS_FACE_CONNECTIVITY_H_
#define BINDINGS_FACE_CONNECTIVITY_H_

#include <pybind11/pybind11.h>

#include "face_connectivity_unstruct.h"

namespace py = pybind11;

namespace bindings
{

/**
 * @brief Bind FaceConnectivityUnstructData and FaceConnectivityUnstruct
 *
 * Should be called BEFORE binding ModelUnstructData (dependency order)
 */
template <typename FloatType, typename ScalarType>
void bindFaceConnectivityUnstruct(py::module& m, const std::string& suffix)
{
  // Bind the Data struct (for Python injection)
  using FaceConnData =
      model::FaceConnectivityUnstructData<FloatType, ScalarType>;

  py::class_<FaceConnData>(m, ("FaceConnectivityUnstructData" + suffix).c_str())
      .def(py::init<>())
      .def_readwrite("n_faces", &FaceConnData::n_faces,
                     "Total number of unique faces")
      .def_readwrite("ndofs_per_face", &FaceConnData::ndofs_per_face,
                     "Number of DOFs (nodes) per face")
      .def_readwrite("elem_to_faces", &FaceConnData::elem_to_faces,
                     "Element to face mapping [n_element x 6]")
      .def_readwrite("face_dofs", &FaceConnData::face_dofs,
                     "Face DOFs [n_faces x ndofs_per_face]")
      .def_readwrite("face_elem_owner", &FaceConnData::face_elem_owner,
                     "Face owner element [n_faces]")
      .def_readwrite("face_elem_neighbor", &FaceConnData::face_elem_neighbor,
                     "Face neighbor element [n_faces], -1 if boundary")
      .def_readwrite("face_local_owner", &FaceConnData::face_local_owner,
                     "Local face index in owner [n_faces]")
      .def_readwrite("face_local_neighbor", &FaceConnData::face_local_neighbor,
                     "Local face index in neighbor [n_faces], -1 if boundary");

  // Bind the class (for runtime usage)
  using FaceConn = model::FaceConnectivityUnstruct<FloatType, ScalarType>;

  py::class_<FaceConn>(m, ("FaceConnectivityUnstruct" + suffix).c_str())
      .def(py::init<>(), "Default constructor")
      .def(py::init<const FaceConnData&>(), py::arg("data"),
           "Construct from data structure")
      .def("build", &FaceConn::build, py::arg("mesh"),
           "Build face connectivity from mesh")
      .def("getNumberOfFaces", &FaceConn::getNumberOfFaces,
           "Get total number of faces")
      .def("getDofsPerFace", &FaceConn::getDofsPerFace,
           "Get number of DOFs per face")
      .def("getGlobalFace", &FaceConn::getGlobalFace, py::arg("elem"),
           py::arg("local_face"),
           "Get global face ID from element and local face")
      .def("getGlobalNodeFromFace", &FaceConn::getGlobalNodeFromFace,
           py::arg("face_id"), py::arg("local_dof"),
           "Get global node index from face and local DOF")
      .def("isBoundaryFace", &FaceConn::isBoundaryFace, py::arg("face_id"),
           "Check if face is on boundary")
      .def("elemOwner", &FaceConn::elemOwner, py::arg("face_id"),
           "Get owner element of face")
      .def("elemNeighbor", &FaceConn::elemNeighbor, py::arg("face_id"),
           "Get neighbor element of face (-1 if boundary)")
      .def("localFaceOwner", &FaceConn::localFaceOwner, py::arg("face_id"),
           "Get local face index in owner element")
      .def("localFaceNeighbor", &FaceConn::localFaceNeighbor,
           py::arg("face_id"), "Get local face index in neighbor element");
}

}  // namespace bindings

#endif  // BINDINGS_FACE_CONNECTIVITY_H_