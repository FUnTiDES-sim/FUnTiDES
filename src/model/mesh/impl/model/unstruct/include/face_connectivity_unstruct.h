#ifndef SRC_MODEL_MESH_MODEL_UNSTRUCTURED_FACE_CONNECTIVITY_UNSTRUCT_H_
#define SRC_MODEL_MESH_MODEL_UNSTRUCTURED_FACE_CONNECTIVITY_UNSTRUCT_H_

#include <model.h>

#include <algorithm>
#include <array>
#include <map>

namespace model
{

/**
 * @brief GPU-compatible face connectivity structure
 */
template <typename ScalarType>
struct FaceConnectivity
{
  ARRAY_INT_VIEW elem_to_faces_;
  ARRAY_INT_VIEW face_dofs_;
  VECTOR_INT_VIEW face_elem_owner_;
  VECTOR_INT_VIEW face_elem_neighbor_;
  VECTOR_INT_VIEW face_local_owner_;
  VECTOR_INT_VIEW face_local_neighbor_;

  int n_faces_ = 0;
  int ndofs_per_face_ = 0;

  PROXY_HOST_DEVICE
  ScalarType numFaces() const { return n_faces_; }

  PROXY_HOST_DEVICE
  ScalarType globalFace(ScalarType elem, CubicFace local_face) const
  {
    return elem_to_faces_(elem, static_cast<int>(local_face));
  }

  PROXY_HOST_DEVICE
  ScalarType globalNode(ScalarType face_id, int local_dof) const
  {
    return face_dofs_(face_id, local_dof);
  }

  PROXY_HOST_DEVICE
  bool isBoundary(ScalarType face_id) const
  {
    return face_elem_neighbor_[face_id] == -1;
  }

  PROXY_HOST_DEVICE
  ScalarType elemOwner(ScalarType face_id) const
  {
    return face_elem_owner_[face_id];
  }

  PROXY_HOST_DEVICE
  ScalarType elemNeighbor(ScalarType face_id) const
  {
    return face_elem_neighbor_[face_id];
  }

  PROXY_HOST_DEVICE
  int localFaceOwner(ScalarType face_id) const
  {
    return face_local_owner_[face_id];
  }

  PROXY_HOST_DEVICE
  int localFaceNeighbor(ScalarType face_id) const
  {
    return face_local_neighbor_[face_id];
  }
};

/**
 * @brief Helper functions for face connectivity construction
 */
namespace FaceConnectivityUtils
{

/**
 * @brief Extract the 4 corner vertices for a face
 */
template <typename ScalarType, typename MeshType>
std::array<ScalarType, 4> extractFaceCorners(const MeshType* mesh,
                                             ScalarType elem,
                                             CubicFace local_face)
{
  switch (local_face)
  {
    case CubicFace::kXMinus:
      return {mesh->globalVertexIndex(elem, 0, 0, 0),
              mesh->globalVertexIndex(elem, 0, 1, 0),
              mesh->globalVertexIndex(elem, 0, 1, 1),
              mesh->globalVertexIndex(elem, 0, 0, 1)};
    case CubicFace::kXPlus:
      return {mesh->globalVertexIndex(elem, 1, 0, 0),
              mesh->globalVertexIndex(elem, 1, 1, 0),
              mesh->globalVertexIndex(elem, 1, 1, 1),
              mesh->globalVertexIndex(elem, 1, 0, 1)};
    case CubicFace::kYMinus:
      return {mesh->globalVertexIndex(elem, 0, 0, 0),
              mesh->globalVertexIndex(elem, 1, 0, 0),
              mesh->globalVertexIndex(elem, 1, 0, 1),
              mesh->globalVertexIndex(elem, 0, 0, 1)};
    case CubicFace::kYPlus:
      return {mesh->globalVertexIndex(elem, 0, 1, 0),
              mesh->globalVertexIndex(elem, 1, 1, 0),
              mesh->globalVertexIndex(elem, 1, 1, 1),
              mesh->globalVertexIndex(elem, 0, 1, 1)};
    case CubicFace::kZMinus:
      return {mesh->globalVertexIndex(elem, 0, 0, 0),
              mesh->globalVertexIndex(elem, 1, 0, 0),
              mesh->globalVertexIndex(elem, 1, 1, 0),
              mesh->globalVertexIndex(elem, 0, 1, 0)};
    case CubicFace::kZPlus:
      return {mesh->globalVertexIndex(elem, 0, 0, 1),
              mesh->globalVertexIndex(elem, 1, 0, 1),
              mesh->globalVertexIndex(elem, 1, 1, 1),
              mesh->globalVertexIndex(elem, 0, 1, 1)};
  }
  return {};
}

/**
 * @brief Create canonical face key (sorted vertices for uniqueness)
 */
template <typename ScalarType>
std::array<ScalarType, 4> makeFaceKey(std::array<ScalarType, 4> corners)
{
  std::sort(corners.begin(), corners.end());
  return corners;
}

}  // namespace FaceConnectivityUtils

}  // namespace model

#endif  // SRC_MODEL_MESH_MODEL_UNSTRUCTURED_FACE_CONNECTIVITY_UNSTRUCT_H_