#ifndef SRC_MODEL_MESH_API_FACE_CONNECTIVITY_H_
#define SRC_MODEL_MESH_API_FACE_CONNECTIVITY_H_

#include "model.h"

namespace model
{

/**
 * @brief Abstract interface for face connectivity
 *
 * Full virtual API — exactement le même pattern que ModelApi.
 * FaceConnectivityStruct implémente tout on-the-fly via formules cartésiennes.
 * FaceConnectivityUnstruct implémente tout via des Kokkos views pré-calculées.
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivityApi
{
 public:
  PROXY_HOST_DEVICE FaceConnectivityApi() = default;
  PROXY_HOST_DEVICE virtual ~FaceConnectivityApi() = default;

  /**
   * @brief Get total number of unique faces in the mesh
   */
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfFaces() const = 0;

  /**
   * @brief Get number of DOFs (nodes) per face
   */
  PROXY_HOST_DEVICE virtual int getDofsPerFace() const = 0;

  /**
   * @brief Get global face ID from element and local face
   * @param elem Element index
   * @param local_face Local face identifier
   * @return Global face ID
   */
  PROXY_HOST_DEVICE virtual ScalarType getGlobalFace(
      ScalarType elem, CubicFace local_face) const = 0;

  /**
   * @brief Get global node index from face and local DOF
   * @param face_id Global face ID
   * @param local_dof Local DOF index on face [0, (order+1)²)
   * @return Global node index
   */
  PROXY_HOST_DEVICE virtual ScalarType getGlobalNodeFromFace(
      ScalarType face_id, int local_dof) const = 0;

  /**
   * @brief Check if face is on domain boundary (no neighbor)
   * @param face_id Global face ID
   * @return True if boundary face
   */
  PROXY_HOST_DEVICE virtual bool isBoundaryFace(ScalarType face_id) const = 0;

  /**
   * @brief Get owner element of a face
   * @param face_id Global face ID
   * @return Owner element index
   */
  PROXY_HOST_DEVICE virtual ScalarType elemOwner(ScalarType face_id) const = 0;

  /**
   * @brief Get neighbor element of a face (-1 if boundary)
   * @param face_id Global face ID
   * @return Neighbor element index, or -1 if boundary
   */
  PROXY_HOST_DEVICE virtual ScalarType elemNeighbor(
      ScalarType face_id) const = 0;

  /**
   * @brief Get local face index of the owner element
   * @param face_id Global face ID
   * @return Local face index (0-5) in owner element
   */
  PROXY_HOST_DEVICE virtual int localFaceOwner(ScalarType face_id) const = 0;

  /**
   * @brief Get local face index of the neighbor element (-1 if boundary)
   * @param face_id Global face ID
   * @return Local face index (0-5) in neighbor element, or -1 if boundary
   */
  PROXY_HOST_DEVICE virtual int localFaceNeighbor(ScalarType face_id) const = 0;
};

}  // namespace model

#endif  // SRC_MODEL_MESH_API_FACE_CONNECTIVITY_H_