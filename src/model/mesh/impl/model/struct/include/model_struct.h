#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_MODEL_STRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_MODEL_STRUCT_H_
#include <elasticity_utils.h>
#include <model.h>

#include "face_connectivity_struct.h"
#include "gllpoints.h"

namespace model
{

/**
 * @brief Data structure for structured Cartesian mesh initialization
 */
template <typename FloatType, typename ScalarType>
struct ModelStructData final : public ModelDataBase<FloatType, ScalarType>
{
 public:
  PROXY_HOST_DEVICE ModelStructData() = default;
  PROXY_HOST_DEVICE ~ModelStructData() = default;
  PROXY_HOST_DEVICE ModelStructData(const ModelStructData&) = default;
  PROXY_HOST_DEVICE ModelStructData& operator=(const ModelStructData&) =
      default;

  ScalarType ex_, ey_, ez_;
  FloatType dx_, dy_, dz_;
  FloatType ox_{0}, oy_{0}, oz_{0};  // Local origin
  VECTOR_REAL_VIEW boundaries_t_;

  bool isModelOnNodes_;
  bool isElastic_;

  /// Optional per-element material arrays. If empty, getters fall back to
  /// hardcoded uniform values. Populated by the builder for heterogeneous
  /// configurations (e.g. acoustoelastic bicouche).
  VECTOR_REAL_VIEW model_vp_element_;   ///< Per-element Vp  (empty → 1500)
  VECTOR_REAL_VIEW model_vs_element_;   ///< Per-element Vs  (empty → 755)
  VECTOR_REAL_VIEW model_rho_element_;  ///< Per-element rho (empty → 1)
  VECTOR_REAL_VIEW model_vp_node_;      ///< Per-node Vp     (empty → 1500)
  VECTOR_REAL_VIEW model_vs_node_;      ///< Per-node Vs     (empty → 755)
  VECTOR_REAL_VIEW model_rho_node_;     ///< Per-node rho    (empty → 1)
};

/**
 * @brief Structured 3D Cartesian mesh with spectral elements
 *
 * Regular hexahedral mesh avec face connectivity entièrement on-the-fly :
 * FaceConnectivityStruct ne stocke que 3 entiers (ex, ey, ez) et calcule
 * tout par formules cartésiennes — zéro allocation supplémentaire.
 *
 * @tparam FloatType Floating point type (float or double)
 * @tparam ScalarType Integer type for indexing
 * @tparam Order Polynomial order of spectral elements
 */
template <typename FloatType, typename ScalarType, int Order>
class ModelStruct : public ModelApi<FloatType, ScalarType>
{
 public:
  using IndexType = std::array<int, 3>;

  PROXY_HOST_DEVICE ModelStruct() = default;

  PROXY_HOST_DEVICE ModelStruct(
      const ModelStructData<FloatType, ScalarType>& data)
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
        free_surface_enabled_(true),
        face_connectivity_(data.ex_, data.ey_, data.ez_, Order),
        model_vp_element_(data.model_vp_element_),
        model_vs_element_(data.model_vs_element_),
        model_rho_element_(data.model_rho_element_),
        model_vp_node_(data.model_vp_node_),
        model_vs_node_(data.model_vs_node_),
        model_rho_node_(data.model_rho_node_)
  {
    nx_ = Order * ex_ + 1;
    ny_ = Order * ey_ + 1;
    nz_ = Order * ez_ + 1;

    hx_ = lx_ / ex_;
    hy_ = ly_ / ey_;
    hz_ = lz_ / ez_;
  }

  PROXY_HOST_DEVICE ModelStruct(const ModelStruct&) = default;
  PROXY_HOST_DEVICE ModelStruct& operator=(const ModelStruct&) = default;
  PROXY_HOST_DEVICE ~ModelStruct() = default;

  // ============================================================================
  // MESH GEOMETRY
  // ============================================================================

  PROXY_HOST_DEVICE
  IndexType elementIndex(const int linearIndex) const
  {
    IndexType elemIndex;
    elemIndex[2] = linearIndex / (ex_ * ey_);
    int const rem = linearIndex - elemIndex[2] * (ex_ * ey_);
    elemIndex[1] = rem / ex_;
    elemIndex[0] = rem - elemIndex[1] * ex_;
    return elemIndex;
  }

  PROXY_HOST_DEVICE
  IndexType globalVertexIndex(IndexType e, int const i, int const j,
                              int const k) const
  {
    return {e[0] + i, e[1] + j, e[2] + k};
  }

  PROXY_HOST_DEVICE
  void vertexCoords(IndexType dofGlobal, FloatType* const coords) const
  {
    coords[0] = dofGlobal[0] * hx_ + ox_;
    coords[1] = dofGlobal[1] * hy_ + oy_;
    coords[2] = dofGlobal[2] * hz_ + oz_;
  }

  PROXY_HOST_DEVICE
  FloatType nodeCoord(ScalarType dofGlobal, int dim) const final
  {
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

    if (localIdx == Order &&
        elemIdx < (dim == 0 ? ex_ : (dim == 1 ? ey_ : ez_)) - 1)
    {
      elemIdx++;
      localIdx = 0;
    }

    FloatType gllPoint = GLLPoints::get(Order, localIdx);
    FloatType elementSize = (dim == 0) ? hx_ : ((dim == 1) ? hy_ : hz_);
    FloatType elementStart = elemIdx * elementSize;
    FloatType physicalCoord =
        elementStart + (gllPoint + 1.0) * elementSize * 0.5;

    switch (dim)
    {
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

  PROXY_HOST_DEVICE
  ScalarType globalNodeIndex(ScalarType e, int i, int j, int k) const final
  {
    ScalarType elemZ = e / (ex_ * ey_);
    ScalarType tmp = e % (ex_ * ey_);
    ScalarType elemY = tmp / ex_;
    ScalarType elemX = tmp % ex_;

    int ix = elemX * Order + i;
    int iy = elemY * Order + j;
    int iz = elemZ * Order + k;

    return ix + iy * nx_ + iz * nx_ * ny_;
  }

  // ============================================================================
  // MATERIAL PROPERTIES
  // ============================================================================

  /// @brief Returns Vp at node @p n. Uses stored array when available, else
  /// 1500.
  PROXY_HOST_DEVICE FloatType getModelVpOnNodes(ScalarType n) const final
  {
    if (model_vp_node_.extent(0) > 0) return model_vp_node_[n];
    return static_cast<FloatType>(1500);
  }
  /// @brief Returns Vp of element @p e. Uses stored array when available, else
  /// 1500.
  PROXY_HOST_DEVICE FloatType getModelVpOnElement(ScalarType e) const final
  {
    if (model_vp_element_.extent(0) > 0) return model_vp_element_[e];
    return static_cast<FloatType>(1500);
  }
  /// @brief Returns rho at node @p n. Uses stored array when available, else 1.
  PROXY_HOST_DEVICE FloatType getModelRhoOnNodes(ScalarType n) const final
  {
    if (model_rho_node_.extent(0) > 0) return model_rho_node_[n];
    return static_cast<FloatType>(1);
  }
  /// @brief Returns rho of element @p e. Uses stored array when available,
  /// else 1.
  PROXY_HOST_DEVICE FloatType getModelRhoOnElement(ScalarType e) const final
  {
    if (model_rho_element_.extent(0) > 0) return model_rho_element_[e];
    return static_cast<FloatType>(1);
  }
  /// @brief Returns Vs at node @p n. Uses stored array when available, else
  /// 755.
  PROXY_HOST_DEVICE FloatType getModelVsOnNodes(ScalarType n) const final
  {
    if (model_vs_node_.extent(0) > 0) return model_vs_node_[n];
    return static_cast<FloatType>(755);
  }
  /// @brief Returns Vs of element @p e. Uses stored array when available, else
  /// 755.
  PROXY_HOST_DEVICE FloatType getModelVsOnElement(ScalarType e) const final
  {
    if (model_vs_element_.extent(0) > 0) return model_vs_element_[e];
    return static_cast<FloatType>(755);
  }
  PROXY_HOST_DEVICE FloatType getModelDeltaOnNodes(ScalarType n) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE FloatType getModelDeltaOnElement(ScalarType e) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnNodes(ScalarType n) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnElement(ScalarType e) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE FloatType getModelGammaOnNodes(ScalarType n) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE FloatType getModelGammaOnElement(ScalarType e) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE ScalarType getModelThetaOnNodes(ScalarType n) const final
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getModelThetaOnElement(ScalarType e) const final
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getModelPhiOnNodes(ScalarType n) const final
  {
    return 0.0;
  }
  PROXY_HOST_DEVICE ScalarType getModelPhiOnElement(ScalarType e) const final
  {
    return 0.0;
  }

  void initElasticityTensors(AnisotropyType anisotropy_type) override
  {
    if (!isElastic_) return;

    if (anisotropy_type == AnisotropyType::kIso ||
        anisotropy_type == AnisotropyType::kVTI)
    {
      // No precomputation needed for ISOTROPIC and VTI, computed on-the-fly
      // inside solver
      return;
    }

    if (anisotropy_type == AnisotropyType::kTTI)
    {
      int n_element = ex_ * ey_ * ez_;

      model_C_tensor_element_ = allocateArray3D<array3DReal>(n_element, 6, 6);

      auto& C_tensor = model_C_tensor_element_;

      MAINLOOPHEAD(n_element, i)
      FloatType CTTI[6][6];

      FloatType vp = 1500.0;
      FloatType vs = 755.0;
      FloatType rho = 1.0;
      FloatType delta = 0.;
      FloatType epsilon = 0.;
      FloatType gamma = 0.0;
      FloatType theta = 0.0;
      FloatType phi = 0.0;

      computeCTensor(vp, vs, rho, delta, epsilon, gamma, theta, phi, CTTI);

      for (int k = 0; k < 6; k++)
        for (int l = 0; l < 6; l++) C_tensor(i, k, l) = CTTI[k][l];
      MAINLOOPEND
    }
  }

  PROXY_HOST_DEVICE
  void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const final
  {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) CTTI[i][j] = model_C_tensor_element_(e, i, j);
  }

  // ============================================================================
  // MESH QUERIES
  // ============================================================================

  PROXY_HOST_DEVICE ScalarType getNumberOfElements() const final
  {
    return ex_ * ey_ * ez_;
  }
  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const final
  {
    return (Order * ex_ + 1) * (Order * ey_ + 1) * (Order * ez_ + 1);
  }
  PROXY_HOST_DEVICE int getNumberOfPointsPerElement() const final
  {
    return (Order + 1) * (Order + 1) * (Order + 1);
  }
  PROXY_HOST_DEVICE int getOrder() const final { return Order; }

  PROXY_HOST_DEVICE
  void faceNormal(ScalarType e, CubicFace local_face,
                  FloatType v[3]) const final
  {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    int direction = static_cast<int>(local_face) / 2;
    FloatType sign = (static_cast<int>(local_face) % 2) ? +1.0 : -1.0;
    v[direction] = sign;
  }

  PROXY_HOST_DEVICE FloatType domainSize(int dim) const final
  {
    switch (dim)
    {
      case 0:
        return lx_;
      case 1:
        return ly_;
      case 2:
        return lz_;
      default:
        return -1;
    }
  }

  PROXY_HOST_DEVICE FloatType getMinSpacing() const final
  {
    if constexpr (Order == 1) return min(hx_, min(hy_, hz_));
    if constexpr (Order == 2) return min(hx_, min(hy_, hz_)) / 2;
    if constexpr (Order == 3) return min(hx_, min(hy_, hz_)) * 0.276393;
    if constexpr (Order == 4) return min(hx_, min(hy_, hz_)) * 0.172673;
    if constexpr (Order == 5) return min(hx_, min(hy_, hz_)) * 0.117472;
    return -1;
  }

  FloatType getMaxSpeed() const final { return 1500; }
  PROXY_HOST_DEVICE bool isModelOnNodes() const final
  {
    return isModelOnNodes_;
  }
  PROXY_HOST_DEVICE bool isElastic() const final { return isElastic_; }

  // ============================================================================
  // BOUNDARY FLAGS
  // ============================================================================

  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(ScalarType n) const override
  {
    if (boundaries_t_.extent(0) == 0) return BoundaryFlag::InteriorNode;
    return static_cast<BoundaryFlag>(boundaries_t_(n));
  }

  PROXY_HOST_DEVICE
  bool isFreeSurface(ScalarType n) const override
  {
    if (boundaries_t_.extent(0) == 0) return false;
    return (boundaries_t_(n) == static_cast<uint8_t>(BoundaryFlag::Surface));
  }

  void initializeBoundaryFlags(bool free_surface_on_top) override
  {
    free_surface_enabled_ = free_surface_on_top;
  }

  void setFreeSurfaceEnabled(bool enable) override
  {
    free_surface_enabled_ = enable;
  }

  void initFreeSurface() override
  {
    if (boundaries_t_.extent(0) > 0) return;  // Déjà pré-calculé par le builder

    boundaries_t_ =
        allocateVector<VECTOR_REAL_VIEW>(getNumberOfNodes(), "boundaries");

    FloatType tol = getMinSpacing() * 1e-4;
    FloatType x_min = ox_, x_max = ox_ + lx_;
    FloatType y_min = oy_, y_max = oy_ + ly_;
    FloatType z_min = oz_, z_max = oz_ + lz_;
    bool enabled = free_surface_enabled_;

    auto boundaries = boundaries_t_;
    auto mesh_copy = *this;

    LOOPHEAD(getNumberOfNodes(), n)
    {
      FloatType x = mesh_copy.nodeCoord(n, 0);
      FloatType y = mesh_copy.nodeCoord(n, 1);
      FloatType z = mesh_copy.nodeCoord(n, 2);

      bool at_xmin = (fabs(x - x_min) < tol);
      bool at_xmax = (fabs(x - x_max) < tol);
      bool at_ymin = (fabs(y - y_min) < tol);
      bool at_ymax = (fabs(y - y_max) < tol);
      bool at_zmin = (fabs(z - z_min) < tol);
      bool at_zmax = (fabs(z - z_max) < tol);

      bool on_boundary =
          at_xmin || at_xmax || at_ymin || at_ymax || at_zmin || at_zmax;

      if (!on_boundary)
        boundaries(n) = static_cast<uint8_t>(BoundaryFlag::InteriorNode);
      else if (at_zmax && enabled)
        boundaries(n) = static_cast<uint8_t>(BoundaryFlag::Surface);
      else
        boundaries(n) = static_cast<uint8_t>(BoundaryFlag::Damping);
    }
    LOOPEND
  }

  /**
   * @brief Build face connectivity
   *
   * Pour ModelStruct : réinitialise FaceConnectivityStruct avec ex/ey/ez.
   * Aucune allocation — tout sera calculé à la volée par formules cartésiennes.
   */
  void buildFaceConnectivity() override
  {
    face_connectivity_ =
        FaceConnectivityStruct<FloatType, ScalarType>(ex_, ey_, ez_, Order);
  }

  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const override
  {
    return face_connectivity_.getNumberOfFaces();
  }

  PROXY_HOST_DEVICE ScalarType
  getGlobalFace(ScalarType elem, CubicFace local_face) const override
  {
    return face_connectivity_.getGlobalFace(elem, local_face);
  }

  PROXY_HOST_DEVICE ScalarType
  getGlobalNodeFromFace(ScalarType face_id, int local_dof) const override
  {
    return face_connectivity_.getGlobalNodeFromFace(face_id, local_dof);
  }

  PROXY_HOST_DEVICE bool isBoundaryFace(ScalarType face_id) const override
  {
    return face_connectivity_.isBoundaryFace(face_id);
  }

  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const
  {
    return face_connectivity_.elemOwner(face_id);
  }

  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const
  {
    return face_connectivity_.elemNeighbor(face_id);
  }

  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const
  {
    return face_connectivity_.localFaceOwner(face_id);
  }

  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const
  {
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
  bool free_surface_enabled_;

  array3DReal model_C_tensor_element_;
  VECTOR_REAL_VIEW boundaries_t_;

  // Optional heterogeneous material arrays (empty → uniform hardcoded values).
  VECTOR_REAL_VIEW model_vp_element_;
  VECTOR_REAL_VIEW model_vs_element_;
  VECTOR_REAL_VIEW model_rho_element_;
  VECTOR_REAL_VIEW model_vp_node_;
  VECTOR_REAL_VIEW model_vs_node_;
  VECTOR_REAL_VIEW model_rho_node_;

  FaceConnectivityStruct<FloatType, ScalarType> face_connectivity_;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_MODEL_STRUCT_H_