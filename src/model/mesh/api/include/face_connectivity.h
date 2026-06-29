#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_FACE_CONNECTIVITY_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_FACE_CONNECTIVITY_H_
#include "model.h"

namespace model {

/**
 * @brief Map a 2D face DOF index to the corresponding 3D element DOF index.
 *
 * The 2D face DOF convention follows the loop order used in
 * FaceConnectivityUnstruct::build(): the fast index is always the first
 * tangential direction and the slow index is the second.
 *
 * @param face   Which of the 6 faces of the hexahedron (CubicFace enum).
 * @param face_dof_2d  2D face DOF index in [0, (order+1)^2).
 * @param order  Polynomial order of the element.
 * @return  3D element-local DOF index in [0, (order+1)^3).
 */
PROXY_HOST_DEVICE int faceLocalToElemLocal(CubicFace face, int face_dof_2d, int order) {
  const int n = order + 1;
  const int u = face_dof_2d % n;
  const int v = face_dof_2d / n;
  switch (face) {
    case CubicFace::kXMinus:
      return u * n + v * n * n;
    case CubicFace::kXPlus:
      return order + u * n + v * n * n;
    case CubicFace::kYMinus:
      return u + v * n * n;
    case CubicFace::kYPlus:
      return u + order * n + v * n * n;
    case CubicFace::kZMinus:
      return u + v * n;
    case CubicFace::kZPlus:
      return u + v * n + order * n * n;
    default:
      return -1;
  }
}

/**
 * @brief Abstract interface for face connectivity
 *
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivityApi {
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
  PROXY_HOST_DEVICE virtual ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const = 0;

  /**
   * @brief Get global node index from face and local DOF
   * @param face_id Global face ID
   * @param local_dof Local DOF index on face [0, (order+1)²)
   * @return Global node index
   */
  PROXY_HOST_DEVICE virtual ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const = 0;

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
  PROXY_HOST_DEVICE virtual ScalarType elemNeighbor(ScalarType face_id) const = 0;

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

  /**
   * @brief Map an owner 2D face DOF index to the corresponding neighbor 2D
   * face DOF index for the same physical node.
   *
   * This permutation is needed by DG flux kernels to match degrees of freedom
   * across an interface when the two adjacent elements index the shared face
   * in different orders.
   *
   * @param face_id   Global face ID.
   * @param owner_dof 2D face DOF index in the owner element [0, (N+1)^2].
   * @return          Corresponding 2D face DOF index in the neighbor element.
   */
  PROXY_HOST_DEVICE virtual int getNeighborFaceDof(ScalarType face_id, int owner_dof) const = 0;

  /**
   * @brief Map an neighbor 2D face DOF index to the corresponding owner 2D
   * face DOF index for the same physical node.
   *
   * This permutation is needed by DG flux kernels to match degrees of freedom
   * across an interface when the two adjacent elements index the shared face
   * in different orders.
   *
   * @param face_id       Global face ID.
   * @param neighbor_dof  2D face DOF index in the neighbor element [0, (N+1)^2].
   * @return              Corresponding 2D face DOF index in the owner element.
   */
  PROXY_HOST_DEVICE virtual int getOwnerFaceDof(ScalarType face_id, int neighbor_dof) const = 0;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_FACE_CONNECTIVITY_H_