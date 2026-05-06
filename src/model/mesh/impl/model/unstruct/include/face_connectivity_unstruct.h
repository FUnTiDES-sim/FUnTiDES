#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
#include <algorithm>
#include <array>
#include <map>

#include "face_connectivity.h"

namespace model {

/**
 * @brief Data structure for unstructured face connectivity initialization
 *
 * Public members allow direct injection from Python (e.g., from HDF5 files).
 * Same pattern as ModelUnstructData.
 */
template <typename FloatType, typename ScalarType>
struct FaceConnectivityUnstructData {
  ScalarType n_faces = 0;
  int ndofs_per_face = 0;

  ARRAY_INT_VIEW elem_to_faces;
  ARRAY_INT_VIEW face_dofs;
  ARRAY_INT_VIEW face_perm;
  VECTOR_INT_VIEW face_elem_owner;
  VECTOR_INT_VIEW face_elem_neighbor;
  VECTOR_INT_VIEW face_local_owner;
  VECTOR_INT_VIEW face_local_neighbor;
};

/**
 * @brief Face connectivity for unstructured meshes
 *
 * Implements FaceConnectivityApi via pre-computed Kokkos views.
 * Can be constructed from FaceConnectivityUnstructData or built from mesh.
 *
 * @tparam FloatType Floating point type
 * @tparam ScalarType Integer type for indexing
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivityUnstruct : public FaceConnectivityApi<FloatType, ScalarType> {
 public:
  FaceConnectivityUnstruct() = default;

  /**
   * @brief Construct from data structure (for Python injection)
   */
  PROXY_HOST_DEVICE
  FaceConnectivityUnstruct(const FaceConnectivityUnstructData<FloatType, ScalarType>& data)
      : n_faces_(data.n_faces),
        ndofs_per_face_(data.ndofs_per_face),
        elem_to_faces_(data.elem_to_faces),
        face_dofs_(data.face_dofs),
        face_perm_(data.face_perm),
        face_elem_owner_(data.face_elem_owner),
        face_elem_neighbor_(data.face_elem_neighbor),
        face_local_owner_(data.face_local_owner),
        face_local_neighbor_(data.face_local_neighbor) {}

  /**
   * @brief Build face connectivity from mesh
   *
   * Extracts faces from elements, identifies unique faces, and fills
   * connectivity tables using a map-based approach.
   */
  void build(const ModelApi<FloatType, ScalarType>& mesh) {
    const ScalarType n_element = mesh.getNumberOfElements();
    const int order = mesh.getOrder();
    const ScalarType max_faces = n_element * 6;
    ndofs_per_face_ = (order + 1) * (order + 1);

    // Temporary arrays at maximum size
    auto elem_to_faces_temp = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    auto face_dofs_temp = allocateArray2D<ARRAY_INT_VIEW>(max_faces, ndofs_per_face_);
    auto face_perm_temp = allocateArray2D<ARRAY_INT_VIEW>(max_faces, ndofs_per_face_);
    auto face_elem_owner_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
    auto face_elem_neighbor_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
    auto face_local_owner_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
    auto face_local_neighbor_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);

    for (ScalarType i = 0; i < max_faces; ++i) face_elem_neighbor_temp(i) = -1;

    using FaceKey = std::array<ScalarType, 4>;
    std::map<FaceKey, ScalarType> face_map;
    ScalarType face_count = 0;

    for (ScalarType elem = 0; elem < n_element; ++elem) {
      for (int lf = 0; lf < 6; ++lf) {
        CubicFace local_face = static_cast<CubicFace>(lf);
        auto corners = extractFaceCorners(mesh, elem, local_face);
        auto face_key = makeFaceKey(corners);

        auto it = face_map.find(face_key);
        if (it == face_map.end()) {
          ScalarType face_id = face_count++;
          face_map[face_key] = face_id;

          fillFaceDofs(mesh, elem, local_face, order,
                       [&](int idx, ScalarType node) { face_dofs_temp(face_id, idx) = node; });

          face_elem_owner_temp(face_id) = elem;
          face_local_owner_temp(face_id) = lf;
          elem_to_faces_temp(elem, lf) = face_id;
        } else {
          ScalarType face_id = it->second;
          face_elem_neighbor_temp(face_id) = elem;
          face_local_neighbor_temp(face_id) = lf;
          elem_to_faces_temp(elem, lf) = face_id;

          // Build permutation: for each owner DOF i, find neighbor DOF j
          // such that both map to the same physical node.
          // ndofs_per_face <= (9+1)^2 = 100 (max order in the codebase).
          constexpr int kMaxDofsPerFace = 100;
          ScalarType neigh_dofs[kMaxDofsPerFace];
          fillFaceDofs(mesh, elem, local_face, order, [&](int idx, ScalarType node) { neigh_dofs[idx] = node; });

          for (int i = 0; i < ndofs_per_face_; ++i) {
            ScalarType owner_node = face_dofs_temp(face_id, i);
            for (int j = 0; j < ndofs_per_face_; ++j) {
              if (neigh_dofs[j] == owner_node) {
                face_perm_temp(face_id, i) = j;
                break;
              }
            }
          }
        }
      }
    }

    // Final allocation at exact size + copy
    n_faces_ = face_count;
    elem_to_faces_ = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    face_dofs_ = allocateArray2D<ARRAY_INT_VIEW>(face_count, ndofs_per_face_);
    face_perm_ = allocateArray2D<ARRAY_INT_VIEW>(face_count, ndofs_per_face_);
    face_elem_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    face_elem_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    face_local_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    face_local_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);

    for (ScalarType elem = 0; elem < n_element; ++elem)
      for (int lf = 0; lf < 6; ++lf) elem_to_faces_(elem, lf) = elem_to_faces_temp(elem, lf);

    for (ScalarType f = 0; f < face_count; ++f) {
      face_elem_owner_(f) = face_elem_owner_temp(f);
      face_elem_neighbor_(f) = face_elem_neighbor_temp(f);
      face_local_owner_(f) = face_local_owner_temp(f);
      face_local_neighbor_(f) = face_local_neighbor_temp(f);
      for (int dof = 0; dof < ndofs_per_face_; ++dof) {
        face_dofs_(f, dof) = face_dofs_temp(f, dof);
        face_perm_(f, dof) = face_perm_temp(f, dof);
      }
    }
  }

  // ==========================================================================
  // FaceConnectivityApi implementation
  // ==========================================================================

  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const override { return n_faces_; }

  PROXY_HOST_DEVICE int getDofsPerFace() const override { return ndofs_per_face_; }

  PROXY_HOST_DEVICE ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const override {
    return elem_to_faces_(elem, static_cast<int>(local_face));
  }

  PROXY_HOST_DEVICE ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const override {
    return face_dofs_(face_id, local_dof);
  }

  PROXY_HOST_DEVICE bool isBoundaryFace(ScalarType face_id) const override {
    return face_elem_neighbor_(face_id) == -1;
  }

  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const override { return face_elem_owner_(face_id); }

  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const override { return face_elem_neighbor_(face_id); }

  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const override { return face_local_owner_(face_id); }

  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const override { return face_local_neighbor_(face_id); }

  PROXY_HOST_DEVICE int getNeighborFaceDof(ScalarType face_id, int owner_dof) const override {
    return face_perm_(face_id, owner_dof);
  }

 private:
  ScalarType n_faces_ = 0;
  int ndofs_per_face_ = 0;

  ARRAY_INT_VIEW elem_to_faces_;
  ARRAY_INT_VIEW face_dofs_;
  ARRAY_INT_VIEW face_perm_;
  VECTOR_INT_VIEW face_elem_owner_;
  VECTOR_INT_VIEW face_elem_neighbor_;
  VECTOR_INT_VIEW face_local_owner_;
  VECTOR_INT_VIEW face_local_neighbor_;

  // Helper methods for build()

  /**
   * @brief Fill face DOFs by iterating over the face nodes and invoking a
   * callback for each (local_idx, global_node) pair.
   */
  template <typename FUNC>
  static void fillFaceDofs(const ModelApi<FloatType, ScalarType>& mesh, ScalarType elem, CubicFace local_face,
                           int order, FUNC&& store) {
    int idx = 0;
    switch (local_face) {
      case CubicFace::kXMinus:
        for (int k = 0; k <= order; ++k)
          for (int j = 0; j <= order; ++j) store(idx++, mesh.globalNodeIndex(elem, 0, j, k));
        break;
      case CubicFace::kXPlus:
        for (int k = 0; k <= order; ++k)
          for (int j = 0; j <= order; ++j) store(idx++, mesh.globalNodeIndex(elem, order, j, k));
        break;
      case CubicFace::kYMinus:
        for (int k = 0; k <= order; ++k)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, i, 0, k));
        break;
      case CubicFace::kYPlus:
        for (int k = 0; k <= order; ++k)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, i, order, k));
        break;
      case CubicFace::kZMinus:
        for (int j = 0; j <= order; ++j)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, i, j, 0));
        break;
      case CubicFace::kZPlus:
        for (int j = 0; j <= order; ++j)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, i, j, order));
        break;
    }
  }

  static std::array<ScalarType, 4> extractFaceCorners(const ModelApi<FloatType, ScalarType>& mesh, ScalarType elem,
                                                      CubicFace local_face) {
    const int o = mesh.getOrder();
    switch (local_face) {
      case CubicFace::kXMinus:
        return {mesh.globalNodeIndex(elem, 0, 0, 0), mesh.globalNodeIndex(elem, 0, o, 0),
                mesh.globalNodeIndex(elem, 0, o, o), mesh.globalNodeIndex(elem, 0, 0, o)};
      case CubicFace::kXPlus:
        return {mesh.globalNodeIndex(elem, o, 0, 0), mesh.globalNodeIndex(elem, o, o, 0),
                mesh.globalNodeIndex(elem, o, o, o), mesh.globalNodeIndex(elem, o, 0, o)};
      case CubicFace::kYMinus:
        return {mesh.globalNodeIndex(elem, 0, 0, 0), mesh.globalNodeIndex(elem, o, 0, 0),
                mesh.globalNodeIndex(elem, o, 0, o), mesh.globalNodeIndex(elem, 0, 0, o)};
      case CubicFace::kYPlus:
        return {mesh.globalNodeIndex(elem, 0, o, 0), mesh.globalNodeIndex(elem, o, o, 0),
                mesh.globalNodeIndex(elem, o, o, o), mesh.globalNodeIndex(elem, 0, o, o)};
      case CubicFace::kZMinus:
        return {mesh.globalNodeIndex(elem, 0, 0, 0), mesh.globalNodeIndex(elem, o, 0, 0),
                mesh.globalNodeIndex(elem, o, o, 0), mesh.globalNodeIndex(elem, 0, o, 0)};
      case CubicFace::kZPlus:
        return {mesh.globalNodeIndex(elem, 0, 0, o), mesh.globalNodeIndex(elem, o, 0, o),
                mesh.globalNodeIndex(elem, o, o, o), mesh.globalNodeIndex(elem, 0, o, o)};
    }
    return {};
  }

  static std::array<ScalarType, 4> makeFaceKey(std::array<ScalarType, 4> corners) {
    std::sort(corners.begin(), corners.end());
    return corners;
  }
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_