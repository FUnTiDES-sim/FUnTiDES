#ifndef SRC_MODEL_MESH_IMPL_MODEL_UNSTRUCT_FACE_CONNECTIVITY_UNSTRUCT_H_
#define SRC_MODEL_MESH_IMPL_MODEL_UNSTRUCT_FACE_CONNECTIVITY_UNSTRUCT_H_

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

  /**
   * @brief Get total number of faces in mesh
   * @return Total number of faces
   */
  PROXY_HOST_DEVICE
  ScalarType numFaces() const { return n_faces_; }

  /**
   * @brief Get global face ID from element and local face
   * @param elem Element index
   * @param local_face Local face identifier (kXMinus, kXPlus, etc.)
   * @return Global face ID
   */
  PROXY_HOST_DEVICE
  ScalarType globalFace(ScalarType elem, CubicFace local_face) const
  {
    return elem_to_faces_(elem, static_cast<int>(local_face));
  }

  /**
   * @brief Get global node index from face and local DOF
   * @param face_id Global face ID
   * @param local_dof Local node index on face [0, ndofs_per_face
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType globalNode(ScalarType face_id, int local_dof) const
  {
    return face_dofs_(face_id, local_dof);
  }

  /**
   * @brief Check if face is on domain boundary
   * @param face_id Global face ID
   * @return True if boundary face
   */
  PROXY_HOST_DEVICE
  bool isBoundary(ScalarType face_id) const
  {
    return face_elem_neighbor_[face_id] == -1;
  }

  /**
   * @brief Get element owner and neighbor for face
   * @param face_id Global face ID
   * @return Element owner and neighbor indices
   */
  PROXY_HOST_DEVICE
  ScalarType elemOwner(ScalarType face_id) const
  {
    return face_elem_owner_[face_id];
  }

  /**
   * @brief Get element neighbor for face
   * @param face_id Global face ID
   * @return Element neighbor index
   */
  PROXY_HOST_DEVICE
  ScalarType elemNeighbor(ScalarType face_id) const
  {
    return face_elem_neighbor_[face_id];
  }

  /**
   * @brief Get local face index in owner element
   * @param face_id Global face ID
   * @return Local face index in owner element
   */
  PROXY_HOST_DEVICE
  int localFaceOwner(ScalarType face_id) const
  {
    return face_local_owner_[face_id];
  }

  /**
   * @brief Get local face index in neighbor element
   * @param face_id Global face ID
   * @return Local face index in neighbor element
   */
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
 * @param mesh Mesh object
 * @param elem Element index
 * @param local_face Local face identifier (kXMinus, kXPlus, etc.)
 * @return Array of 4 corner vertex indices
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

#endif  // SRC_MODEL_MESH_IMPL_MODEL_UNSTRUCT_FACE_CONNECTIVITY_UNSTRUCT_H_