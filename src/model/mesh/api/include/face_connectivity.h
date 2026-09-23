#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_FACE_CONNECTIVITY_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_FACE_CONNECTIVITY_H_
#include "model.h"

namespace model {

/**
 * @brief Map a 2D face DOF index to the corresponding element-local DOF index.
 *
 * @see docs/design.md, section "Hexahedron local numbering", for both index conventions.
 * @param[in] face Face of the hexahedron.
 * @param[in] face_dof_2d 2D face DOF index, in [0, (order+1)^2).
 * @param[in] order Polynomial order of the element.
 * @return Element-local DOF index in [0, (order+1)^3), or -1 for an invalid face.
 */
PROXY_HOST_DEVICE constexpr int faceLocalToElemLocal(CubicFace face, int face_dof_2d, int order) {
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
 * @brief Map a 2D face DOF index and a depth along the face normal to the corresponding
 * element-local DOF index.
 *
 * depth is the local index along the normal axis of the face, counted from the minus side:
 * 0 lies on the minus face and order on the plus face, whichever of the two faces is given.
 * With depth equal to the face's own index (0 for minus faces, order for plus faces) the result
 * equals faceLocalToElemLocal(face, face_dof_2d, order).
 *
 * @see docs/design.md, section "Hexahedron local numbering".
 * @param[in] face Face of the hexahedron.
 * @param[in] face_dof_2d 2D face DOF index, in [0, (order+1)^2).
 * @param[in] depth Local index along the face normal axis, in [0, order].
 * @param[in] order Polynomial order of the element.
 * @return Element-local DOF index in [0, (order+1)^3), or -1 for an invalid face.
 */
PROXY_HOST_DEVICE constexpr int faceLocalToElemLocalAtDepth(CubicFace face, int face_dof_2d, int depth, int order) {
  const int n = order + 1;
  const int u = face_dof_2d % n;
  const int v = face_dof_2d / n;
  switch (face) {
    case CubicFace::kXMinus:
    case CubicFace::kXPlus:
      return depth + u * n + v * n * n;
    case CubicFace::kYMinus:
    case CubicFace::kYPlus:
      return u + depth * n + v * n * n;
    case CubicFace::kZMinus:
    case CubicFace::kZPlus:
      return u + v * n + depth * n * n;
    default:
      return -1;
  }
}

/**
 * @brief Abstract face-based view of a hexahedral mesh, used by DG flux kernels.
 *
 * Every geometric face has one global id, one owner element and at most one neighbor element;
 * which of the two adjacent elements is the owner is implementation-defined. Face DOFs follow
 * the owner element's face numbering; getNeighborFaceDof() and getOwnerFaceDof() translate
 * between the two sides.
 * @see docs/design.md, sections "Hexahedron local numbering" and "Device calls on mesh objects".
 *
 * @tparam FloatType Floating-point type, unused by the interface itself.
 * @tparam ScalarType Integer type of element, node and face indices.
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivityApi {
 public:
  PROXY_HOST_DEVICE FaceConnectivityApi() = default;
  PROXY_HOST_DEVICE virtual ~FaceConnectivityApi() = default;

  /**
   * @brief Get the number of distinct faces of the mesh.
   */
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfFaces() const = 0;

  /**
   * @brief Get the number of nodes per face, (order+1)^2.
   */
  PROXY_HOST_DEVICE virtual int getDofsPerFace() const = 0;

  /**
   * @brief Get the global face id of a local face of an element.
   * @param[in] elem Element index.
   * @param[in] local_face Local face of elem.
   * @return Global face id, in [0, getNumberOfFaces()).
   */
  PROXY_HOST_DEVICE virtual ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const = 0;

  /**
   * @brief Get the global node index of a node of a face.
   * @param[in] face_id Global face id.
   * @param[in] local_dof 2D face DOF index in the owner element numbering, in
   * [0, getDofsPerFace()).
   * @return Global node index.
   */
  PROXY_HOST_DEVICE virtual ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const = 0;

  /**
   * @brief Tell whether a face has no neighbor element.
   * @param[in] face_id Global face id.
   */
  PROXY_HOST_DEVICE virtual bool isBoundaryFace(ScalarType face_id) const = 0;

  /**
   * @brief Get the owner element of a face.
   * @param[in] face_id Global face id.
   * @return Owner element index.
   */
  PROXY_HOST_DEVICE virtual ScalarType elemOwner(ScalarType face_id) const = 0;

  /**
   * @brief Get the neighbor element of a face.
   * @param[in] face_id Global face id.
   * @return Neighbor element index, or -1 for a boundary face.
   */
  PROXY_HOST_DEVICE virtual ScalarType elemNeighbor(ScalarType face_id) const = 0;

  /**
   * @brief Get the local face of the owner element that coincides with a face.
   * @param[in] face_id Global face id.
   * @return CubicFace value, in [0, 5].
   */
  PROXY_HOST_DEVICE virtual int localFaceOwner(ScalarType face_id) const = 0;

  /**
   * @brief Get the local face of the neighbor element that coincides with a face.
   * @param[in] face_id Global face id.
   * @return CubicFace value in [0, 5], or -1 for a boundary face.
   */
  PROXY_HOST_DEVICE virtual int localFaceNeighbor(ScalarType face_id) const = 0;

  /**
   * @brief Map a 2D face DOF index of the owner element to the neighbor's 2D face DOF index of
   * the same physical node.
   *
   * Needed because the two adjacent elements may number the shared face in different orders.
   * Meaningful only for faces that have a neighbor.
   *
   * @param[in] face_id Global face id.
   * @param[in] owner_dof 2D face DOF index in the owner element, in [0, getDofsPerFace()).
   * @return 2D face DOF index in the neighbor element.
   */
  PROXY_HOST_DEVICE virtual int getNeighborFaceDof(ScalarType face_id, int owner_dof) const = 0;

  /**
   * @brief Inverse of getNeighborFaceDof(): map a 2D face DOF index of the neighbor element to
   * the owner's 2D face DOF index of the same physical node.
   *
   * @param[in] face_id Global face id.
   * @param[in] neighbor_dof 2D face DOF index in the neighbor element, in [0, getDofsPerFace()).
   * @return 2D face DOF index in the owner element.
   */
  PROXY_HOST_DEVICE virtual int getOwnerFaceDof(ScalarType face_id, int neighbor_dof) const = 0;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_FACE_CONNECTIVITY_H_
