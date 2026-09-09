#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
#include <Kokkos_UnorderedMap.hpp>
#include <limits>
#include <stdexcept>

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

  arrayInt elem_to_faces;
  arrayInt face_dofs;
  arrayInt face_perm;
  arrayInt face_perm_inv;
  vectorInt face_elem_owner;
  vectorInt face_elem_neighbor;
  vectorInt face_local_owner;
  vectorInt face_local_neighbor;
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
template <typename FloatType, typename ScalarType, int ORDER = -1>
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
        face_perm_inv_(data.face_perm_inv),
        face_elem_owner_(data.face_elem_owner),
        face_elem_neighbor_(data.face_elem_neighbor),
        face_local_owner_(data.face_local_owner),
        face_local_neighbor_(data.face_local_neighbor) {}

  /**
   * @brief Canonical (sorted) 4-corner key identifying a face, independent
   * of which adjacent element/local-face it's seen from.
   *
   * Public: nvcc requires types captured by an extended __device__ lambda
   * (used inside build()) to have public accessibility.
   */
  struct FaceKey {
    ScalarType nodes[4];

    KOKKOS_INLINE_FUNCTION bool operator==(const FaceKey& other) const {
      return nodes[0] == other.nodes[0] && nodes[1] == other.nodes[1] && nodes[2] == other.nodes[2] &&
             nodes[3] == other.nodes[3];
    }
  };

  /**
   * @brief Build face connectivity from mesh
   *
   * Extracts faces from elements, identifies unique faces, and fills
   * connectivity tables using a thread-safe map-based approach.
   *
   * Runs as Kokkos::parallel_for kernels (both ModelStruct and ModelUnstruct
   * expose a device-callable globalNodeIndex()) instead of a single-threaded
   * host loop, since the per-face work (map insertion + O(ndofs^2)
   * permutation search) is what dominates build() cost, not memory
   * locality. Face ownership is elected deterministically (element with the
   * smaller index owns the face) via an atomic minimum over a packed
   * (elem, local_face) code, matching the outcome of the original serial
   * elem-ascending loop.
   *
   * @tparam MESH_TYPE Concrete mesh type (struct or unstruct); must expose
   *   a device-callable globalNodeIndex(elem,i,j,k), getNumberOfElements(),
   *   getOrder().
   */
  template <typename MESH_TYPE>
  void build(const MESH_TYPE& mesh) {
    const ScalarType n_element = mesh.getNumberOfElements();
    const int mesh_order = mesh.getOrder();
    const int order = (ORDER >= 0) ? ORDER : mesh_order;
    const ScalarType max_faces = n_element * 6;
    ndofs_per_face_ = (order + 1) * (order + 1);

    using FaceMap = Kokkos::UnorderedMap<FaceKey, void>;
    FaceMap face_map(static_cast<uint32_t>(max_faces));

    // owner_code/face_id_of_bucket/face_count_dev are vectorInt (int-valued
    // Kokkos views, independent of ScalarType) — use int's own max as the
    // "unset" sentinel, not ScalarType's (which may be wider, e.g. long, and
    // would silently truncate through deep_copy into an incorrect value).
    vectorInt owner_code = allocateVector<vectorInt>(face_map.capacity());
    Kokkos::deep_copy(owner_code, std::numeric_limits<int>::max());

    // Pass A: insert the (sorted-corner) key for every element face and
    // atomically elect the owner as the smaller of the (at most two)
    // elements touching it, via a packed (elem, local_face) code.
    Kokkos::parallel_for(
        "FaceConnectivityUnstruct_insert", n_element, KOKKOS_LAMBDA(const ScalarType elem) {
          for (int lf = 0; lf < 6; ++lf) {
            const FaceKey key = makeFaceKey(mesh, mesh_order, elem, static_cast<CubicFace>(lf));
            const auto res = face_map.insert(key);
            Kokkos::atomic_fetch_min(&owner_code(res.index()), static_cast<int>(elem * 8 + lf));
          }
        });
    Kokkos::fence();
    if (face_map.failed_insert()) {
      throw std::runtime_error("FaceConnectivityUnstruct::build: face map insertion failed (capacity too small)");
    }

    // Pass B: compact the sparse map slots into dense face ids [0, face_count).
    // The prefix sum makes ids a pure function of the mesh; an atomic counter
    // would order them by thread scheduling, so two instances built from the
    // same mesh would disagree — and callers share face id lists across
    // instances (see DGsolver::setFaceConnectivity).
    vectorInt face_id_of_bucket = allocateVector<vectorInt>(face_map.capacity());
    ScalarType face_count = 0;
    Kokkos::parallel_scan(
        "FaceConnectivityUnstruct_compact", n_element * 6,
        KOKKOS_LAMBDA(const ScalarType flat, ScalarType& running_id, const bool is_final) {
          const ScalarType elem = flat / 6;
          const int lf = static_cast<int>(flat % 6);
          const FaceKey key = makeFaceKey(mesh, mesh_order, elem, static_cast<CubicFace>(lf));
          const uint32_t idx = face_map.find(key);
          if (owner_code(idx) != static_cast<int>(elem * 8 + lf)) return;
          if (is_final) face_id_of_bucket(idx) = running_id;
          ++running_id;
        },
        face_count);
    Kokkos::fence();

    // Final device allocation at exact size.
    n_faces_ = face_count;
    elem_to_faces_ = allocateArray2D<arrayInt>(n_element, 6);
    face_dofs_ = allocateArray2D<arrayInt>(face_count, ndofs_per_face_);
    face_perm_ = allocateArray2D<arrayInt>(face_count, ndofs_per_face_);
    face_perm_inv_ = allocateArray2D<arrayInt>(face_count, ndofs_per_face_);
    face_elem_owner_ = allocateVector<vectorInt>(face_count);
    face_elem_neighbor_ = allocateVector<vectorInt>(face_count);
    face_local_owner_ = allocateVector<vectorInt>(face_count);
    face_local_neighbor_ = allocateVector<vectorInt>(face_count);
    Kokkos::deep_copy(face_elem_neighbor_, -1);

    arrayInt face_dofs = face_dofs_;
    vectorInt face_elem_owner = face_elem_owner_;
    vectorInt face_local_owner = face_local_owner_;
    arrayInt elem_to_faces = elem_to_faces_;
    const int ndofs_per_face = ndofs_per_face_;

    // Pass C: record elem->face for every element; the elected owner fills
    // face_dofs_ and owner metadata.
    Kokkos::parallel_for(
        "FaceConnectivityUnstruct_owner", n_element, KOKKOS_LAMBDA(const ScalarType elem) {
          for (int lf = 0; lf < 6; ++lf) {
            const FaceKey key = makeFaceKey(mesh, mesh_order, elem, static_cast<CubicFace>(lf));
            const uint32_t idx = face_map.find(key);
            const ScalarType face_id = face_id_of_bucket(idx);
            elem_to_faces(elem, lf) = face_id;
            if (owner_code(idx) == static_cast<int>(elem * 8 + lf)) {
              fillFaceDofs(mesh, elem, static_cast<CubicFace>(lf), order,
                           [&](int d, ScalarType node) { face_dofs(face_id, d) = node; });
              face_elem_owner(face_id) = elem;
              face_local_owner(face_id) = lf;
            }
          }
        });
    Kokkos::fence();

    arrayInt face_perm = face_perm_;
    arrayInt face_perm_inv = face_perm_inv_;
    vectorInt face_elem_neighbor = face_elem_neighbor_;
    vectorInt face_local_neighbor = face_local_neighbor_;

    // Pass D: the non-owner side fills neighbor metadata and the
    // owner<->neighbor DOF permutation (reads face_dofs_ written in Pass C).
    Kokkos::parallel_for(
        "FaceConnectivityUnstruct_neighbor", n_element, KOKKOS_LAMBDA(const ScalarType elem) {
          for (int lf = 0; lf < 6; ++lf) {
            const ScalarType face_id = elem_to_faces(elem, lf);
            if (face_elem_owner(face_id) == elem && face_local_owner(face_id) == lf) continue;

            face_elem_neighbor(face_id) = elem;
            face_local_neighbor(face_id) = lf;

            // ndofs_per_face <= (9+1)^2 = 100 (max order in the codebase).
            constexpr int kMaxDofsPerFace = 100;
            ScalarType neigh_dofs[kMaxDofsPerFace];
            fillFaceDofs(mesh, elem, static_cast<CubicFace>(lf), order,
                         [&](int d, ScalarType node) { neigh_dofs[d] = node; });

            for (int i = 0; i < ndofs_per_face; ++i) {
              const ScalarType owner_node = face_dofs(face_id, i);
              for (int j = 0; j < ndofs_per_face; ++j) {
                if (neigh_dofs[j] == owner_node) {
                  face_perm(face_id, i) = j;
                  face_perm_inv(face_id, j) = i;
                  break;
                }
              }
            }
          }
        });
    // Callers read the connectivity from the host right after build()
    // (list construction, tagging); Pass D must be complete before they do.
    Kokkos::fence();
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

  PROXY_HOST_DEVICE int getOwnerFaceDof(ScalarType face_id, int neighbor_dof) const override {
    return face_perm_inv_(face_id, neighbor_dof);
  }

 private:
  ScalarType n_faces_ = 0;
  int ndofs_per_face_ = 0;

  arrayInt elem_to_faces_;
  arrayInt face_dofs_;
  arrayInt face_perm_;
  arrayInt face_perm_inv_;
  vectorInt face_elem_owner_;
  vectorInt face_elem_neighbor_;
  vectorInt face_local_owner_;
  vectorInt face_local_neighbor_;

  // Helper methods for build()

  /**
   * @brief Fill face DOFs by iterating over the face nodes and invoking a
   * callback for each (local_idx, global_node) pair.
   */
  template <typename MESH_TYPE, typename FUNC>
  KOKKOS_INLINE_FUNCTION static void fillFaceDofs(const MESH_TYPE& mesh, ScalarType elem, CubicFace local_face,
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

  /**
   * @brief Build the 4-corner key identifying a local face.
   *
   * Extracts the 4 corner global node indices of the given local face, always
   * in the same per-orientation traversal order. Two elements sharing a
   * physical face see it from opposite orientations (e.g. kXPlus/kXMinus)
   * whose traversal orders line up on the same global corner nodes, so both
   * produce an equal FaceKey without needing to sort the 4 nodes. This
   * relies on globalNodeIndex() returning identical global indices for
   * shared corners regardless of which adjacent element queries them, and
   * is what makes the Kokkos::UnorderedMap-based matching in build() work.
   *
   * @tparam MESH_TYPE Mesh type exposing a device-callable globalNodeIndex().
   * @param mesh Mesh providing node indexing.
   * @param mesh_order Polynomial order of the mesh.
   * @param elem Element index owning the local face.
   * @param local_face Local face identifier (see CubicFace).
   * @return FaceKey holding the 4 corner global node indices of the face.
   */
  template <typename MESH_TYPE>
  KOKKOS_INLINE_FUNCTION static FaceKey makeFaceKey(const MESH_TYPE& mesh, int mesh_order, ScalarType elem,
                                                    CubicFace local_face) {
    const int o = mesh_order;
    FaceKey key{};
    switch (local_face) {
      case CubicFace::kXMinus:
        key = {mesh.globalNodeIndex(elem, 0, 0, 0), mesh.globalNodeIndex(elem, 0, o, 0),
               mesh.globalNodeIndex(elem, 0, o, o), mesh.globalNodeIndex(elem, 0, 0, o)};
        break;
      case CubicFace::kXPlus:
        key = {mesh.globalNodeIndex(elem, o, 0, 0), mesh.globalNodeIndex(elem, o, o, 0),
               mesh.globalNodeIndex(elem, o, o, o), mesh.globalNodeIndex(elem, o, 0, o)};
        break;
      case CubicFace::kYMinus:
        key = {mesh.globalNodeIndex(elem, 0, 0, 0), mesh.globalNodeIndex(elem, o, 0, 0),
               mesh.globalNodeIndex(elem, o, 0, o), mesh.globalNodeIndex(elem, 0, 0, o)};
        break;
      case CubicFace::kYPlus:
        key = {mesh.globalNodeIndex(elem, 0, o, 0), mesh.globalNodeIndex(elem, o, o, 0),
               mesh.globalNodeIndex(elem, o, o, o), mesh.globalNodeIndex(elem, 0, o, o)};
        break;
      case CubicFace::kZMinus:
        key = {mesh.globalNodeIndex(elem, 0, 0, 0), mesh.globalNodeIndex(elem, o, 0, 0),
               mesh.globalNodeIndex(elem, o, o, 0), mesh.globalNodeIndex(elem, 0, o, 0)};
        break;
      case CubicFace::kZPlus:
        key = {mesh.globalNodeIndex(elem, 0, 0, o), mesh.globalNodeIndex(elem, o, 0, o),
               mesh.globalNodeIndex(elem, o, o, o), mesh.globalNodeIndex(elem, 0, o, o)};
        break;
    }
    // Sort the 4 corners so both adjacent elements produce the same key
    // regardless of local orientation (tiny fixed-size network, unrolled).
    ScalarType* n = key.nodes;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3 - a; ++b)
        if (n[b] > n[b + 1]) {
          const ScalarType tmp = n[b];
          n[b] = n[b + 1];
          n[b + 1] = tmp;
        }
    return key;
  }
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
