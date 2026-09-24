#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_FACE_CONNECTIVITY_UNSTRUCT_H_
#include <Kokkos_UnorderedMap.hpp>
#include <limits>
#include <stdexcept>

#include "face_connectivity.h"

namespace model {

/**
 * @brief Plain data used to initialize a FaceConnectivityUnstruct from pre-computed tables.
 *
 * Members are public so the tables can be injected directly (for example from Python).
 * Table layouts are those of the matching members of FaceConnectivityUnstruct.
 */
template <typename FloatType, typename ScalarType>
struct FaceConnectivityUnstructData {
  ScalarType n_faces = 0;         ///< Number of unique faces.
  int ndofs_per_face = 0;         ///< Nodes per face, (order + 1)^2.
  arrayInt elem_to_faces;         ///< Global face id, shape (numElements, 6), indexed by CubicFace.
  arrayInt face_dofs;             ///< Global node index of each owner-side face dof, shape (n_faces, ndofs_per_face).
  arrayInt face_perm;             ///< Owner dof to neighbor dof, shape (n_faces, ndofs_per_face).
  arrayInt face_perm_inv;         ///< Neighbor dof to owner dof, shape (n_faces, ndofs_per_face).
  vectorInt face_elem_owner;      ///< Owner element, size n_faces.
  vectorInt face_elem_neighbor;   ///< Neighbor element, size n_faces, -1 on a boundary face.
  vectorInt face_local_owner;     ///< Local face index seen from the owner, size n_faces.
  vectorInt face_local_neighbor;  ///< Local face index seen from the neighbor, size n_faces.
};

/**
 * @brief Face connectivity of an unstructured hexahedral mesh, stored in Kokkos views.
 *
 * Can be filled from a FaceConnectivityUnstructData or built from a mesh with build().
 * Every face is shared by at most two elements: the one with the smaller index is the
 * owner, the other one the neighbor.
 *
 * @tparam FloatType Floating point type.
 * @tparam ScalarType Integer type used for element, face and node indices.
 * @tparam ORDER Polynomial order of the face dofs; -1 means the order of the mesh given to build().
 */
template <typename FloatType, typename ScalarType, int ORDER = -1>
class FaceConnectivityUnstruct : public FaceConnectivityApi<FloatType, ScalarType> {
 public:
  FaceConnectivityUnstruct() = default;

  /**
   * @brief Construct from pre-computed tables.
   * @param data Tables to copy (the views are shared, not deep-copied).
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
   * @brief Key identifying a face by its 4 corner nodes, sorted in ascending order.
   *
   * Independent of the element and local face the face is seen from.
   * Public because nvcc requires types captured by an extended __device__ lambda
   * (used in build()) to be publicly accessible.
   */
  struct FaceKey {
    ScalarType nodes[4];  ///< Global node indices of the 4 corners, ascending.

    KOKKOS_INLINE_FUNCTION bool operator==(const FaceKey& other) const {
      return nodes[0] == other.nodes[0] && nodes[1] == other.nodes[1] && nodes[2] == other.nodes[2] &&
             nodes[3] == other.nodes[3];
    }
  };

  /**
   * @brief Build the face connectivity of a mesh.
   *
   * Identifies the unique faces and fills all tables on the device. Face ids are a
   * deterministic function of the mesh (ascending element, then local face, of the owner).
   * The tables are complete and visible to the host when the call returns.
   *
   * @tparam MESH_TYPE Mesh type; must provide a device-callable globalNodeIndex(elem, i, j, k),
   *   getNumberOfElements() and getOrder().
   * @param[in] mesh Mesh to analyze.
   * @param[in] geom_order Geometric order of the node grid of the mesh, used to locate the "Plus"
   *   faces and to rescale the face dof indices. Pass it only when ORDER is lower than the order
   *   of the mesh (p-adaptive coupling); -1 samples the face dofs at ORDER.
   * @throws std::runtime_error If the face map overflows.
   */
  template <typename MESH_TYPE>
  void build(const MESH_TYPE& mesh, int geom_order = -1) {
    const ScalarType n_element = mesh.getNumberOfElements();
    const int mesh_order = mesh.getOrder();
    const int order = (ORDER >= 0) ? ORDER : mesh_order;
    const ScalarType max_faces = n_element * 6;
    ndofs_per_face_ = (order + 1) * (order + 1);

    using FaceMap = Kokkos::UnorderedMap<FaceKey, void>;
    FaceMap face_map(static_cast<uint32_t>(max_faces));

    // These views are vectorInt, independent of ScalarType: the "unset" sentinel is
    // int's max, because a wider ScalarType would be truncated by deep_copy.
    vectorInt owner_code = allocateVector<vectorInt>(face_map.capacity());
    Kokkos::deep_copy(owner_code, std::numeric_limits<int>::max());

    // Pass A: insert every element face and elect the owner as the smaller of the
    // (at most two) packed (elem, local_face) codes touching it.
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
    // A prefix sum is used instead of an atomic counter so that ids do not depend on
    // thread scheduling: callers share face id lists across instances built from the
    // same mesh (see DGsolver::setFaceConnectivity).
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

    // Pass C: every element records its face ids; the owner also fills face_dofs_
    // and the owner metadata.
    Kokkos::parallel_for(
        "FaceConnectivityUnstruct_owner", n_element, KOKKOS_LAMBDA(const ScalarType elem) {
          for (int lf = 0; lf < 6; ++lf) {
            const FaceKey key = makeFaceKey(mesh, mesh_order, elem, static_cast<CubicFace>(lf));
            const uint32_t idx = face_map.find(key);
            const ScalarType face_id = face_id_of_bucket(idx);
            elem_to_faces(elem, lf) = face_id;
            if (owner_code(idx) == static_cast<int>(elem * 8 + lf)) {
              fillFaceDofs(
                  mesh, elem, static_cast<CubicFace>(lf), order,
                  [&](int d, ScalarType node) { face_dofs(face_id, d) = node; }, geom_order);
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

    // Pass D: the non-owner side fills the neighbor metadata and the dof permutations.
    // Needs face_dofs_ from Pass C.
    Kokkos::parallel_for(
        "FaceConnectivityUnstruct_neighbor", n_element, KOKKOS_LAMBDA(const ScalarType elem) {
          for (int lf = 0; lf < 6; ++lf) {
            const ScalarType face_id = elem_to_faces(elem, lf);
            if (face_elem_owner(face_id) == elem && face_local_owner(face_id) == lf) continue;

            face_elem_neighbor(face_id) = elem;
            face_local_neighbor(face_id) = lf;

            // Upper bound: (9 + 1)^2 dofs at the maximum order of the codebase.
            constexpr int kMaxDofsPerFace = 100;
            ScalarType neigh_dofs[kMaxDofsPerFace];
            fillFaceDofs(
                mesh, elem, static_cast<CubicFace>(lf), order, [&](int d, ScalarType node) { neigh_dofs[d] = node; },
                geom_order);

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
    // Callers read the tables from the host right after build().
    Kokkos::fence();
  }

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
  ScalarType n_faces_ = 0;  ///< Number of unique faces.
  int ndofs_per_face_ = 0;  ///< Nodes per face, (order + 1)^2.

  arrayInt elem_to_faces_;         ///< Global face id, shape (numElements, 6), indexed by CubicFace.
  arrayInt face_dofs_;             ///< Global node of each owner-side face dof, shape (n_faces, ndofs_per_face).
  arrayInt face_perm_;             ///< Owner dof to neighbor dof, shape (n_faces, ndofs_per_face).
  arrayInt face_perm_inv_;         ///< Neighbor dof to owner dof, shape (n_faces, ndofs_per_face).
  vectorInt face_elem_owner_;      ///< Owner element, size n_faces.
  vectorInt face_elem_neighbor_;   ///< Neighbor element, size n_faces, -1 on a boundary face.
  vectorInt face_local_owner_;     ///< Local face index seen from the owner, size n_faces.
  vectorInt face_local_neighbor_;  ///< Local face index seen from the neighbor, size n_faces.

  /**
   * @brief Enumerate the nodes of a local face and pass each to a callback.
   *
   * @tparam MESH_TYPE Mesh type providing a device-callable globalNodeIndex(elem, i, j, k).
   * @tparam FUNC Callable with signature void(int local_dof, ScalarType global_node).
   * @param[in] mesh Mesh providing node indexing.
   * @param[in] elem Element index.
   * @param[in] local_face Local face of the element.
   * @param[in] order Polynomial order of the face dofs; ndofs_per_face_ = (order + 1)^2.
   * @param[in] store Callback invoked once per face dof, in storage order.
   * @param[in] geom_order Geometric order of the node grid of the element. Used for the fixed
   *   face-normal coordinate of the "Plus" faces and to rescale the tangential indices onto that
   *   grid, so that a lower-order sub-solver sharing a higher-order mesh still sees the real
   *   element boundary. -1 means @p order.
   */
  template <typename MESH_TYPE, typename FUNC>
  KOKKOS_INLINE_FUNCTION static void fillFaceDofs(const MESH_TYPE& mesh, ScalarType elem, CubicFace local_face,
                                                  int order, FUNC&& store, int geom_order = -1) {
    int const far = (geom_order >= 0) ? geom_order : order;
    auto rescale = [&](int i) { return (geom_order >= 0 && order > 0) ? (i * geom_order) / order : i; };
    int idx = 0;
    switch (local_face) {
      case CubicFace::kXMinus:
        for (int k = 0; k <= order; ++k)
          for (int j = 0; j <= order; ++j) store(idx++, mesh.globalNodeIndex(elem, 0, rescale(j), rescale(k)));
        break;
      case CubicFace::kXPlus:
        for (int k = 0; k <= order; ++k)
          for (int j = 0; j <= order; ++j) store(idx++, mesh.globalNodeIndex(elem, far, rescale(j), rescale(k)));
        break;
      case CubicFace::kYMinus:
        for (int k = 0; k <= order; ++k)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, rescale(i), 0, rescale(k)));
        break;
      case CubicFace::kYPlus:
        for (int k = 0; k <= order; ++k)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, rescale(i), far, rescale(k)));
        break;
      case CubicFace::kZMinus:
        for (int j = 0; j <= order; ++j)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, rescale(i), rescale(j), 0));
        break;
      case CubicFace::kZPlus:
        for (int j = 0; j <= order; ++j)
          for (int i = 0; i <= order; ++i) store(idx++, mesh.globalNodeIndex(elem, rescale(i), rescale(j), far));
        break;
    }
  }

  /**
   * @brief Build the key of a local face from its 4 corner nodes.
   *
   * Two elements sharing a face produce equal keys because globalNodeIndex() returns the
   * same global index for a shared corner whichever element queries it, and the corners
   * are sorted before returning.
   *
   * @tparam MESH_TYPE Mesh type providing a device-callable globalNodeIndex().
   * @param[in] mesh Mesh providing node indexing.
   * @param[in] mesh_order Polynomial order of the mesh.
   * @param[in] elem Element owning the local face.
   * @param[in] local_face Local face of the element.
   * @return Key holding the 4 sorted corner node indices.
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
    // Sort the 4 corners so that the key does not depend on the local orientation.
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
