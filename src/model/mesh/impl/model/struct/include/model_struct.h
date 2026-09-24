#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_MODEL_STRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_MODEL_STRUCT_H_
#include <elasticity_utils.h>
#include <model.h>

#include "face_connectivity_struct.h"
#include "gllpoints.h"

namespace model {

/**
 * @brief Parameters and optional material arrays needed to build a ModelStruct.
 *
 * Filled by the builder of a structured Cartesian mesh, then copied into the model.
 *
 * @tparam FloatType Floating point type.
 * @tparam ScalarType Integer type used for indexing.
 */
template <typename FloatType, typename ScalarType>
struct ModelStructData final : public ModelDataBase<FloatType, ScalarType> {
 public:
  PROXY_HOST_DEVICE ModelStructData() = default;
  PROXY_HOST_DEVICE ~ModelStructData() = default;
  PROXY_HOST_DEVICE ModelStructData(const ModelStructData&) = default;
  PROXY_HOST_DEVICE ModelStructData& operator=(const ModelStructData&) = default;

  ScalarType ex_, ey_, ez_;  ///< Number of elements along x, y, z in the local subdomain.
  /// @todo VERIFY: ModelStruct stores these as the local domain lengths (lx_, ly_, lz_) and derives
  /// the element size from them; is dx_, dy_, dz_ a domain length or an element size?
  FloatType dx_, dy_, dz_;
  FloatType ox_{0.0f}, oy_{0.0f}, oz_{0.0f};  ///< Origin of the local subdomain.
  vectorInt boundaries_t_;  ///< Per-node BoundaryFlag values; empty when no boundary data is provided.

  bool isModelOnNodes_;  ///< True if material properties are given per node, false if per element.
  bool isElastic_;       ///< True for an elastic model, false for an acoustic one.

  /// Optional material arrays. When an array is empty, the matching getter of ModelStruct
  /// returns a hardcoded uniform value (given in each line).
  vectorReal model_vp_element_;   ///< Per-element Vp (empty: 1500).
  vectorReal model_vs_element_;   ///< Per-element Vs (empty: 755).
  vectorReal model_rho_element_;  ///< Per-element rho (empty: 1).
  vectorReal model_vp_node_;      ///< Per-node Vp (empty: 1500).
  vectorReal model_vs_node_;      ///< Per-node Vs (empty: 755).
  vectorReal model_rho_node_;     ///< Per-node rho (empty: 1).
  vectorReal model_qp_element_;   ///< Per-element P quality factor (empty: 1.0e9).
  vectorReal model_qs_element_;   ///< Per-element S quality factor (empty: 1.0e9).
  vectorReal model_qp_node_;      ///< Per-node P quality factor (empty: 1.0e9).
  vectorReal model_qs_node_;      ///< Per-node S quality factor (empty: 1.0e9).
};

/**
 * @brief Structured Cartesian hexahedral mesh with spectral elements.
 *
 * Node and element numbering follow x fastest, then y, then z. All geometry and face
 * connectivity are computed from formulas on the element counts, so the mesh stores no
 * coordinate or connectivity table.
 *
 * @tparam FloatType Floating point type.
 * @tparam ScalarType Integer type used for indexing.
 * @tparam Order Polynomial order of the spectral elements.
 */
template <typename FloatType, typename ScalarType, int Order>
class ModelStruct final : public ModelApi<FloatType, ScalarType> {
 public:
  using IndexType = std::array<int, 3>;

  PROXY_HOST_DEVICE ModelStruct() = default;

  PROXY_HOST_DEVICE ModelStruct(const ModelStructData<FloatType, ScalarType>& data)
      : ex_(data.ex_),
        ey_(data.ey_),
        ez_(data.ez_),
        ox_(data.ox_),
        oy_(data.oy_),
        oz_(data.oz_),
        lx_(data.dx_),
        ly_(data.dy_),
        lz_(data.dz_),
        isModelOnNodes_(data.isModelOnNodes_),
        boundaries_t_(data.boundaries_t_),
        isElastic_(data.isElastic_),
        face_connectivity_(data.ex_, data.ey_, data.ez_, Order),
        model_vp_element_(data.model_vp_element_),
        model_vs_element_(data.model_vs_element_),
        model_rho_element_(data.model_rho_element_),
        model_vp_node_(data.model_vp_node_),
        model_vs_node_(data.model_vs_node_),
        model_rho_node_(data.model_rho_node_),
        model_qp_element_(data.model_qp_element_),
        model_qs_element_(data.model_qs_element_),
        model_qp_node_(data.model_qp_node_),
        model_qs_node_(data.model_qs_node_) {
    nx_ = Order * ex_ + 1;
    ny_ = Order * ey_ + 1;
    nz_ = Order * ez_ + 1;

    hx_ = lx_ / static_cast<FloatType>(ex_);
    hy_ = ly_ / static_cast<FloatType>(ey_);
    hz_ = lz_ / static_cast<FloatType>(ez_);
  }

  PROXY_HOST_DEVICE ModelStruct(const ModelStruct&) = default;
  PROXY_HOST_DEVICE ModelStruct& operator=(const ModelStruct&) = default;
  PROXY_HOST_DEVICE ~ModelStruct() = default;

  /**
   * @brief Converts a linear element index into its (ix, iy, iz) position in the element grid.
   * @param linearIndex Element index, equal to ix + iy * ex + iz * ex * ey.
   * @return Array {ix, iy, iz}.
   */
  PROXY_HOST_DEVICE
  IndexType elementIndex(const int linearIndex) const {
    IndexType elemIndex;
    elemIndex[2] = linearIndex / (ex_ * ey_);
    int const rem = linearIndex - elemIndex[2] * (ex_ * ey_);
    elemIndex[1] = rem / ex_;
    elemIndex[0] = rem - elemIndex[1] * ex_;
    return elemIndex;
  }

  /**
   * @brief Position in the vertex grid of a vertex of an element.
   * @param e Element position {ix, iy, iz} in the element grid.
   * @param i,j,k Vertex offset along x, y, z inside the element (0 or 1).
   * @return Array {ix + i, iy + j, iz + k}.
   */
  PROXY_HOST_DEVICE
  IndexType globalVertexIndex(IndexType e, int const i, int const j, int const k) const {
    return {e[0] + i, e[1] + j, e[2] + k};
  }

  /**
   * @brief Physical coordinates of a vertex of the element grid, origin included.
   * @param dofGlobal Vertex position {ix, iy, iz} in the vertex grid.
   * @param[out] coords Array of 3 values x, y, z.
   */
  PROXY_HOST_DEVICE
  void vertexCoords(IndexType dofGlobal, FloatType* const coords) const {
    coords[0] = static_cast<FloatType>(dofGlobal[0]) * hx_ + ox_;
    coords[1] = static_cast<FloatType>(dofGlobal[1]) * hy_ + oy_;
    coords[2] = static_cast<FloatType>(dofGlobal[2]) * hz_ + oz_;
  }

  /**
   * @brief Coordinate of a global node along one axis, at its GLL position, origin included.
   * @param dofGlobal Global node index, equal to i + j * nx + k * nx * ny.
   * @param dim Axis: 0 for x, 1 for y, 2 for z.
   */
  PROXY_HOST_DEVICE
  FloatType nodeCoord(ScalarType dofGlobal, int dim) const final {
    int nodesPerDim[3];
    nodesPerDim[0] = (ex_ * Order) + 1;
    nodesPerDim[1] = (ey_ * Order) + 1;
    nodesPerDim[2] = (ez_ * Order) + 1;

    int k = dofGlobal / (nodesPerDim[0] * nodesPerDim[1]);
    int remainder = dofGlobal % (nodesPerDim[0] * nodesPerDim[1]);
    int j = remainder / nodesPerDim[0];
    int i = remainder % nodesPerDim[0];

    int nodeIdx[3] = {i, j, k};
    int elemIdx = nodeIdx[dim] / Order;
    int localIdx = nodeIdx[dim] % Order;

    if (localIdx == Order && elemIdx < (dim == 0 ? ex_ : (dim == 1 ? ey_ : ez_)) - 1) {
      elemIdx++;
      localIdx = 0;
    }

    FloatType gllPoint = static_cast<FloatType>(GLLPoints::get(Order, localIdx));
    FloatType elementSize = (dim == 0) ? hx_ : ((dim == 1) ? hy_ : hz_);
    FloatType elementStart = static_cast<FloatType>(elemIdx) * elementSize;
    FloatType physicalCoord = elementStart + (gllPoint + 1.0f) * elementSize * 0.5f;

    switch (dim) {
      case 0:
        physicalCoord += ox_;
        break;
      case 1:
        physicalCoord += oy_;
        break;
      case 2:
        physicalCoord += oz_;
        break;
    }
    return physicalCoord;
  }

  /**
   * @brief Global node index of a local node of an element.
   * @param e Linear element index.
   * @param i,j,k Local node position along x, y, z, each in [0, Order].
   * @return Global node index i_g + j_g * nx + k_g * nx * ny.
   */
  PROXY_HOST_DEVICE
  ScalarType globalNodeIndex(ScalarType e, int i, int j, int k) const final {
    ScalarType elemZ = e / (ex_ * ey_);
    ScalarType tmp = e % (ex_ * ey_);
    ScalarType elemY = tmp / ex_;
    ScalarType elemX = tmp % ex_;
    int ix = elemX * Order + i;
    int iy = elemY * Order + j;
    int iz = elemZ * Order + k;
    return ix + iy * nx_ + iz * nx_ * ny_;
  }

  /// @brief Vp at node @p n; 1500 when no per-node array is stored.
  PROXY_HOST_DEVICE FloatType getModelVpOnNodes(ScalarType n) const final {
    if (model_vp_node_.extent(0) > 0) return model_vp_node_[n];
    return 1500.0f;
  }
  /// @brief Vp of element @p e; 1500 when no per-element array is stored.
  PROXY_HOST_DEVICE FloatType getModelVpOnElement(ScalarType e) const final {
    if (model_vp_element_.extent(0) > 0) return model_vp_element_[e];
    return 1500.0f;
  }
  /// @brief Density at node @p n; 1 when no per-node array is stored.
  PROXY_HOST_DEVICE FloatType getModelRhoOnNodes(ScalarType n) const final {
    if (model_rho_node_.extent(0) > 0) return model_rho_node_[n];
    return 1.0f;
  }
  /// @brief Density of element @p e; 1 when no per-element array is stored.
  PROXY_HOST_DEVICE FloatType getModelRhoOnElement(ScalarType e) const final {
    if (model_rho_element_.extent(0) > 0) return model_rho_element_[e];
    return 1.0f;
  }
  /// @brief Vs at node @p n; 755 when no per-node array is stored.
  PROXY_HOST_DEVICE FloatType getModelVsOnNodes(ScalarType n) const final {
    if (model_vs_node_.extent(0) > 0) return model_vs_node_[n];
    return 755.0f;
  }
  /// @brief Vs of element @p e; 755 when no per-element array is stored.
  PROXY_HOST_DEVICE FloatType getModelVsOnElement(ScalarType e) const final {
    if (model_vs_element_.extent(0) > 0) return model_vs_element_[e];
    return 755.0f;
  }

  /**
   * @brief Overwrites the stored material properties at one node.
   *
   * Only the arrays that are allocated are written. Host-only.
   *
   * @param n Global node index.
   * @param vp P-wave velocity (m/s).
   * @param vs S-wave velocity (m/s).
   * @param rho Density (kg/m^3).
   */
  void setModelNodeProps(ScalarType n, FloatType vp, FloatType vs, FloatType rho) {
    if (model_vp_node_.extent(0) > 0) model_vp_node_[n] = vp;
    if (model_vs_node_.extent(0) > 0) model_vs_node_[n] = vs;
    if (model_rho_node_.extent(0) > 0) model_rho_node_[n] = rho;
  }

  /// @brief P quality factor at node @p n; 1.0e9 when no per-node array is stored.
  PROXY_HOST_DEVICE FloatType getModelQpOnNodes(ScalarType n) const final {
    if (model_qp_node_.extent(0) > 0) return model_qp_node_[n];
    return 1.0e9f;
  }
  /// @brief P quality factor of element @p e; 1.0e9 when no per-element array is stored.
  PROXY_HOST_DEVICE FloatType getModelQpOnElement(ScalarType e) const final {
    if (model_qp_element_.extent(0) > 0) return model_qp_element_[e];
    return 1.0e9f;
  }
  /// @brief S quality factor at node @p n; 1.0e9 when no per-node array is stored.
  PROXY_HOST_DEVICE FloatType getModelQsOnNodes(ScalarType n) const final {
    if (model_qs_node_.extent(0) > 0) return model_qs_node_[n];
    return 1.0e9f;
  }
  /// @brief S quality factor of element @p e; 1.0e9 when no per-element array is stored.
  PROXY_HOST_DEVICE FloatType getModelQsOnElement(ScalarType e) const final {
    if (model_qs_element_.extent(0) > 0) return model_qs_element_[e];
    return 1.0e9f;
  }

  /// @brief Anisotropy parameters: this model has no storage for them and always returns 0.
  PROXY_HOST_DEVICE FloatType getModelDeltaOnNodes(ScalarType n) const final { return 0.0f; }
  PROXY_HOST_DEVICE FloatType getModelDeltaOnElement(ScalarType e) const final { return 0.0f; }
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnNodes(ScalarType n) const final { return 0.0f; }
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnElement(ScalarType e) const final { return 0.0f; }
  PROXY_HOST_DEVICE FloatType getModelGammaOnNodes(ScalarType n) const final { return 0.0f; }
  PROXY_HOST_DEVICE FloatType getModelGammaOnElement(ScalarType e) const final { return 0.0f; }
  PROXY_HOST_DEVICE ScalarType getModelThetaOnNodes(ScalarType n) const final { return 0; }
  PROXY_HOST_DEVICE ScalarType getModelThetaOnElement(ScalarType e) const final { return 0; }
  PROXY_HOST_DEVICE ScalarType getModelPhiOnNodes(ScalarType n) const final { return 0; }
  PROXY_HOST_DEVICE ScalarType getModelPhiOnElement(ScalarType e) const final { return 0; }

  /**
   * @brief Allocates and fills the per-element 6x6 stiffness tensors.
   *
   * Does nothing for an acoustic model or for isotropic and VTI anisotropy. For TTI, allocates
   * one 6x6 tensor per element.
   */
  void initElasticityTensors(AnisotropyType anisotropy_type) override {
    if (!isElastic_) return;
    if (anisotropy_type == AnisotropyType::kIso || anisotropy_type == AnisotropyType::kVTI) return;

    if (anisotropy_type == AnisotropyType::kTTI) {
      int n_element = ex_ * ey_ * ez_;
      model_C_tensor_element_ = allocateArray3D<array3DReal>(n_element, 6, 6);
      auto C_tensor = model_C_tensor_element_;

      Kokkos::parallel_for(
          "Model init ElasticTensors",
          Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_element),
          KOKKOS_LAMBDA(const int i) {
            FloatType CTTI[6][6];
            FloatType vp = 1500.0f, vs = 755.0f, rho = 1.0f;
            FloatType delta = 0.0f, epsilon = 0.0f, gamma = 0.0f;
            FloatType theta = 0.0f, phi = 0.0f;
            computeCTensor(vp, vs, rho, delta, epsilon, gamma, theta, phi, CTTI);
            for (int k = 0; k < 6; k++)
              for (int l = 0; l < 6; l++) C_tensor(i, k, l) = CTTI[k][l];
          });
    }
  }

  /**
   * @brief Copies the stored 6x6 stiffness tensor of an element.
   * @pre initElasticityTensors() allocated the tensors (TTI only).
   * @param e Element index.
   * @param[out] CTTI 6x6 tensor.
   */
  PROXY_HOST_DEVICE void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const final {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) CTTI[i][j] = model_C_tensor_element_(e, i, j);
  }

  PROXY_HOST_DEVICE ScalarType getNumberOfElements() const final { return ex_ * ey_ * ez_; }
  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const final { return nx_ * ny_ * nz_; }
  PROXY_HOST_DEVICE int getNumberOfPointsPerElement() const final { return (Order + 1) * (Order + 1) * (Order + 1); }
  PROXY_HOST_DEVICE int getOrder() const final { return Order; }

  /**
   * @brief Outward unit normal of a face of an element.
   *
   * The face index is assumed to encode axis = local_face / 2 and sign = -1 for even, +1 for odd.
   *
   * @param e Element index (unused: the mesh is axis-aligned).
   * @param local_face Local face of the element.
   * @param[out] v Array of 3 components.
   */
  PROXY_HOST_DEVICE void faceNormal(ScalarType e, CubicFace local_face, FloatType v[3]) const final {
    v[0] = 0.0f;
    v[1] = 0.0f;
    v[2] = 0.0f;
    int direction = static_cast<int>(local_face) / 2;
    FloatType sign = (static_cast<int>(local_face) % 2) ? 1.0f : -1.0f;
    v[direction] = sign;
  }

  /// @brief Length of the local subdomain along axis @p dim (0, 1, 2), or -1 for another value.
  PROXY_HOST_DEVICE FloatType domainSize(int dim) const final {
    if (dim == 0) return lx_;
    if (dim == 1) return ly_;
    if (dim == 2) return lz_;
    return -1.0f;
  }

  /**
   * @brief Smallest distance between two neighbouring GLL nodes.
   *
   * Element size times the smallest normalized GLL spacing of the order. Returns -1 for an
   * order above 9.
   */
  PROXY_HOST_DEVICE FloatType getMinSpacing() const final {
    const FloatType h = min(hx_, min(hy_, hz_));
    if constexpr (Order == 1) return h;
    if constexpr (Order == 2) return h * 0.5000000000f;
    if constexpr (Order == 3) return h * 0.2763932023f;
    if constexpr (Order == 4) return h * 0.1726731646f;
    if constexpr (Order == 5) return h * 0.1174723380f;
    if constexpr (Order == 6) return h * 0.0848880519f;
    if constexpr (Order == 7) return h * 0.0641299257f;
    if constexpr (Order == 8) return h * 0.0501210023f;
    if constexpr (Order == 9) return h * 0.0402330459f;
    return -1.0f;
  }

  FloatType getMaxSpeed() const final { return 1500.0f; }
  PROXY_HOST_DEVICE bool isModelOnNodes() const final { return isModelOnNodes_; }
  PROXY_HOST_DEVICE bool isElastic() const final { return isElastic_; }

  /// @brief Boundary flag of node @p n; InteriorNode when no boundary data is stored.
  PROXY_HOST_DEVICE BoundaryFlag boundaryType(ScalarType n) const override {
    if (boundaries_t_.extent(0) == 0) return BoundaryFlag::InteriorNode;
    return static_cast<BoundaryFlag>(boundaries_t_(n));
  }

  /// @brief True if node @p n is flagged Surface; false when no boundary data is stored.
  PROXY_HOST_DEVICE bool isFreeSurface(ScalarType n) const override {
    if (boundaries_t_.extent(0) == 0) return false;
    return boundaries_t_(n) == static_cast<ScalarType>(BoundaryFlag::Surface);
  }

  /// @brief Sets a uniform per-element quality factor for P and S waves; host-only.
  void setQualityFactors(FloatType qp, FloatType qs) override {
    ScalarType nElem = getNumberOfElements();
    model_qp_element_ = allocateVector<vectorReal>(nElem, "model_qp_element");
    model_qs_element_ = allocateVector<vectorReal>(nElem, "model_qs_element");
    for (ScalarType e = 0; e < nElem; ++e) {
      model_qp_element_(e) = qp;
      model_qs_element_(e) = qs;
    }
  }

  /// @brief Rebuilds the face connectivity from the element counts.
  void buildFaceConnectivity() override {
    face_connectivity_ = FaceConnectivityStruct<FloatType, ScalarType>(ex_, ey_, ez_, Order);
  }

  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const override { return face_connectivity_.getNumberOfFaces(); }
  PROXY_HOST_DEVICE ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const override {
    return face_connectivity_.getGlobalFace(elem, local_face);
  }
  PROXY_HOST_DEVICE ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const override {
    return face_connectivity_.getGlobalNodeFromFace(face_id, local_dof);
  }

  /**
   * @brief True if the face lies on the domain boundary.
   *
   * With boundary data, a face is a boundary face when none of its nodes is an InteriorNode.
   * Without it, the answer comes from the face connectivity.
   */
  PROXY_HOST_DEVICE bool isBoundaryFace(ScalarType face_id) const override {
    if (boundaries_t_.extent(0) == 0) return face_connectivity_.isBoundaryFace(face_id);
    int const n_dofs = face_connectivity_.getDofsPerFace();
    for (int q = 0; q < n_dofs; ++q) {
      if (boundaries_t_(getGlobalNodeFromFace(face_id, q)) == static_cast<ScalarType>(BoundaryFlag::InteriorNode))
        return false;
    }
    return true;
  }

  /// @brief Face queries forwarded to the face connectivity.
  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const { return face_connectivity_.elemOwner(face_id); }
  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const {
    return face_connectivity_.elemNeighbor(face_id);
  }
  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const { return face_connectivity_.localFaceOwner(face_id); }
  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const {
    return face_connectivity_.localFaceNeighbor(face_id);
  }

 private:
  ScalarType ex_, ey_, ez_;
  ScalarType nx_, ny_, nz_;
  FloatType lx_, ly_, lz_;
  FloatType hx_, hy_, hz_;
  FloatType ox_, oy_, oz_;
  bool isModelOnNodes_;
  bool isElastic_;
  array3DReal model_C_tensor_element_;
  vectorInt boundaries_t_;
  vectorReal model_qp_element_, model_qs_element_;
  vectorReal model_qp_node_, model_qs_node_;
  vectorReal model_vp_element_, model_vs_element_, model_rho_element_;
  vectorReal model_vp_node_, model_vs_node_, model_rho_node_;
  FaceConnectivityStruct<FloatType, ScalarType> face_connectivity_;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_MODEL_STRUCT_H_
