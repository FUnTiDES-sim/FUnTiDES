#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_FACE_CONNECTIVITY_STRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_FACE_CONNECTIVITY_STRUCT_H_
#include "face_connectivity.h"

namespace model {

/**
 * @brief Face connectivity of a structured Cartesian hexahedral mesh, computed from
 * the element counts alone (no tables are stored).
 *
 * Elements are numbered with i fastest, then j, then k. Global faces are numbered
 * X faces first, then Y faces, then Z faces. Global nodes are numbered
 * ix + iy * nx + iz * nx * ny, with nx = order * ex + 1 and ny = order * ey + 1.
 *
 * @tparam FloatType Floating point type
 * @tparam ScalarType Integer type for indexing
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivityStruct : public FaceConnectivityApi<FloatType, ScalarType> {
 public:
  PROXY_HOST_DEVICE FaceConnectivityStruct() = default;

  /**
   * @brief Build the connectivity of an ex x ey x ez element mesh.
   * @param[in] ex Number of elements along x
   * @param[in] ey Number of elements along y
   * @param[in] ez Number of elements along z
   * @param[in] order Polynomial order of the elements
   */
  PROXY_HOST_DEVICE
  FaceConnectivityStruct(ScalarType ex, ScalarType ey, ScalarType ez, int order)
      : ex_(ex), ey_(ey), ez_(ez), order_(order) {
    nx_ = order_ * ex_ + 1;
    ny_ = order_ * ey_ + 1;
    offset_y_ = (ex_ + 1) * ey_ * ez_;
    offset_z_ = offset_y_ + ex_ * (ey_ + 1) * ez_;
  }

  /** @brief Total number of faces (X, Y and Z faces). */
  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const override { return offset_z_ + ex_ * ey_ * (ez_ + 1); }

  /** @brief Number of nodes on one face, (order + 1)^2. */
  PROXY_HOST_DEVICE int getDofsPerFace() const override { return (order_ + 1) * (order_ + 1); }

  /**
   * @brief Global id of a face of an element.
   * @param[in] elem Element index
   * @param[in] local_face Local face of the element
   * @return Global face id, or -1 if local_face is not a valid face
   */
  PROXY_HOST_DEVICE ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const override {
    ScalarType elem_k = elem / (ex_ * ey_);
    ScalarType tmp = elem % (ex_ * ey_);
    ScalarType elem_j = tmp / ex_;
    ScalarType elem_i = tmp % ex_;

    switch (local_face) {
      case CubicFace::kXMinus:
        return elem_i + elem_j * (ex_ + 1) + elem_k * (ex_ + 1) * ey_;
      case CubicFace::kXPlus:
        return (elem_i + 1) + elem_j * (ex_ + 1) + elem_k * (ex_ + 1) * ey_;
      case CubicFace::kYMinus:
        return offset_y_ + elem_i + elem_j * ex_ + elem_k * ex_ * (ey_ + 1);
      case CubicFace::kYPlus:
        return offset_y_ + elem_i + (elem_j + 1) * ex_ + elem_k * ex_ * (ey_ + 1);
      case CubicFace::kZMinus:
        return offset_z_ + elem_i + elem_j * ex_ + elem_k * ex_ * ey_;
      case CubicFace::kZPlus:
        return offset_z_ + elem_i + elem_j * ex_ + (elem_k + 1) * ex_ * ey_;
      default:
        return -1;
    }
  }

  /**
   * @brief Global node index of a node of a face.
   *
   * The local dof runs over the two tangential directions, the first one fastest:
   * - X face: local_dof = k * (order + 1) + j
   * - Y face: local_dof = k * (order + 1) + i
   * - Z face: local_dof = j * (order + 1) + i
   *
   * @param[in] face_id Global face id
   * @param[in] local_dof Local dof on the face, in [0, (order + 1)^2)
   * @return Global node index
   */
  PROXY_HOST_DEVICE ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const override {
    ScalarType ix, iy, iz;

    if (face_id < offset_y_)  // X face
    {
      ScalarType i_face = face_id % (ex_ + 1);
      ScalarType j_face = (face_id / (ex_ + 1)) % ey_;
      ScalarType k_face = face_id / ((ex_ + 1) * ey_);

      ScalarType j_local = local_dof % (order_ + 1);
      ScalarType k_local = local_dof / (order_ + 1);

      ix = i_face * order_;
      iy = j_face * order_ + j_local;
      iz = k_face * order_ + k_local;
    } else if (face_id < offset_z_)  // Y face
    {
      ScalarType local = face_id - offset_y_;
      ScalarType i_face = local % ex_;
      ScalarType j_face = (local / ex_) % (ey_ + 1);
      ScalarType k_face = local / (ex_ * (ey_ + 1));

      ScalarType i_local = local_dof % (order_ + 1);
      ScalarType k_local = local_dof / (order_ + 1);

      ix = i_face * order_ + i_local;
      iy = j_face * order_;
      iz = k_face * order_ + k_local;
    } else  // Z face
    {
      ScalarType local = face_id - offset_z_;
      ScalarType i_face = local % ex_;
      ScalarType j_face = (local / ex_) % ey_;
      ScalarType k_face = local / (ex_ * ey_);

      ScalarType i_local = local_dof % (order_ + 1);
      ScalarType j_local = local_dof / (order_ + 1);

      ix = i_face * order_ + i_local;
      iy = j_face * order_ + j_local;
      iz = k_face * order_;
    }

    return ix + iy * nx_ + iz * nx_ * ny_;
  }

  /**
   * @brief Tell whether a face lies on the mesh boundary (has no neighbor).
   * @param[in] face_id Global face id
   * @return True if the face is on the outer boundary of the mesh
   */
  PROXY_HOST_DEVICE bool isBoundaryFace(ScalarType face_id) const override {
    if (face_id < offset_y_) {
      ScalarType i = face_id % (ex_ + 1);
      return (i == 0 || i == ex_);
    } else if (face_id < offset_z_) {
      ScalarType local = face_id - offset_y_;
      ScalarType j = (local / ex_) % (ey_ + 1);
      return (j == 0 || j == ey_);
    } else {
      ScalarType local = face_id - offset_z_;
      ScalarType k = local / (ex_ * ey_);
      return (k == 0 || k == ez_);
    }
  }

  /**
   * @brief Owner element of a face: the element on its minus side, or the only
   * adjacent element for a face on the maximum boundary.
   * @param[in] face_id Global face id
   * @return Owner element index
   */
  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const override {
    if (face_id < offset_y_) {
      ScalarType i = face_id % (ex_ + 1);
      ScalarType j = (face_id / (ex_ + 1)) % ey_;
      ScalarType k = face_id / ((ex_ + 1) * ey_);
      ScalarType ei = (i < ex_) ? i : i - 1;
      return ei + j * ex_ + k * ex_ * ey_;
    } else if (face_id < offset_z_) {
      ScalarType local = face_id - offset_y_;
      ScalarType i = local % ex_;
      ScalarType j = (local / ex_) % (ey_ + 1);
      ScalarType k = local / (ex_ * (ey_ + 1));
      ScalarType ej = (j < ey_) ? j : j - 1;
      return i + ej * ex_ + k * ex_ * ey_;
    } else {
      ScalarType local = face_id - offset_z_;
      ScalarType i = local % ex_;
      ScalarType j = (local / ex_) % ey_;
      ScalarType k = local / (ex_ * ey_);
      ScalarType ek = (k < ez_) ? k : k - 1;
      return i + j * ex_ + ek * ex_ * ey_;
    }
  }

  /**
   * @brief Neighbor element of a face, the one on the opposite side from the owner.
   * @param[in] face_id Global face id
   * @return Neighbor element index, or -1 if the face is on the boundary
   */
  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const override {
    if (isBoundaryFace(face_id)) return -1;

    if (face_id < offset_y_) {
      ScalarType i = face_id % (ex_ + 1);
      ScalarType j = (face_id / (ex_ + 1)) % ey_;
      ScalarType k = face_id / ((ex_ + 1) * ey_);
      return (i - 1) + j * ex_ + k * ex_ * ey_;
    } else if (face_id < offset_z_) {
      ScalarType local = face_id - offset_y_;
      ScalarType i = local % ex_;
      ScalarType j = (local / ex_) % (ey_ + 1);
      ScalarType k = local / (ex_ * (ey_ + 1));
      return i + (j - 1) * ex_ + k * ex_ * ey_;
    } else {
      ScalarType local = face_id - offset_z_;
      ScalarType i = local % ex_;
      ScalarType j = (local / ex_) % ey_;
      ScalarType k = local / (ex_ * ey_);
      return i + j * ex_ + (k - 1) * ex_ * ey_;
    }
  }

  /**
   * @brief Local face index of a face in its owner element.
   * @param[in] face_id Global face id
   * @return Value of the CubicFace enumerator, in [0, 5]
   */
  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const override {
    if (face_id < offset_y_) {
      ScalarType i = face_id % (ex_ + 1);
      return (i < ex_) ? static_cast<int>(CubicFace::kXMinus) : static_cast<int>(CubicFace::kXPlus);
    } else if (face_id < offset_z_) {
      ScalarType local = face_id - offset_y_;
      ScalarType j = (local / ex_) % (ey_ + 1);
      return (j < ey_) ? static_cast<int>(CubicFace::kYMinus) : static_cast<int>(CubicFace::kYPlus);
    } else {
      ScalarType local = face_id - offset_z_;
      ScalarType k = local / (ex_ * ey_);
      return (k < ez_) ? static_cast<int>(CubicFace::kZMinus) : static_cast<int>(CubicFace::kZPlus);
    }
  }

  /**
   * @brief Local face index of a face in its neighbor element.
   * @param[in] face_id Global face id
   * @return Value of the CubicFace enumerator, or -1 if the face is on the boundary
   */
  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const override {
    if (isBoundaryFace(face_id)) return -1;
    return localFaceOwner(face_id) ^ 1;  // opposite face; relies on the CubicFace enumerator order
  }

  /**
   * @brief Face dof of the neighbor matching a face dof of the owner.
   *
   * Adjacent elements index a shared face in the same order, so this is the identity.
   */
  PROXY_HOST_DEVICE int getNeighborFaceDof(ScalarType /*face_id*/, int owner_dof) const override { return owner_dof; }

  /**
   * @brief Face dof of the owner matching a face dof of the neighbor.
   *
   * Adjacent elements index a shared face in the same order, so this is the identity.
   */
  PROXY_HOST_DEVICE int getOwnerFaceDof(ScalarType /*face_id*/, int neighbor_dof) const override {
    return neighbor_dof;
  }

 private:
  ScalarType ex_{0}, ey_{0}, ez_{0};      ///< Number of elements along x, y, z
  ScalarType nx_{0}, ny_{0};              ///< Number of nodes along x and y
  ScalarType offset_y_{0}, offset_z_{0};  ///< Id of the first Y face and of the first Z face
  int order_{0};                          ///< Polynomial order of the elements
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_FACE_CONNECTIVITY_STRUCT_H_
