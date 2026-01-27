#ifndef SRC_MODEL_MESH_IMPL_MODEL_STRUCT_MODEL_STRUCT_H_
#define SRC_MODEL_MESH_IMPL_MODEL_STRUCT_MODEL_STRUCT_H_

#include <elasticity_utils.h>
#include <model.h>

#include "gllpoints.h"

namespace model
{

/**
 * @brief Data structure for structured Cartesian mesh initialization
 *
 * Contains grid dimensions and domain size for constructing a regular mesh.
 *
 * @tparam FloatType Floating point type for coordinates and properties
 * @tparam ScalarType Integer type for indices and counts
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
  FloatType ox_{0}, oy_{0}, oz_{0};

  bool isModelOnNodes_;
  bool isElastic_;
};

/**
 * @brief Structured 3D Cartesian mesh with spectral elements
 *
 * Regular hexahedral mesh with implicit connectivity computed via formulas.
 * Optimized for uniform domains with O(1) face lookup and zero memory overhead
 * for face connectivity.
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

  /**
   * @brief Default constructor
   */
  PROXY_HOST_DEVICE ModelStruct() = default;

  /**
   * @brief Construct from data structure
   * @param data Mesh configuration parameters
   */
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
        isElastic_(data.isElastic_),
        freeSurfaceEnabled_(true)
  {
    nx_ = Order * ex_ + 1;
    ny_ = Order * ey_ + 1;
    nz_ = Order * ez_ + 1;

    hx_ = lx_ / ex_;
    hy_ = ly_ / ey_;
    hz_ = lz_ / ez_;

    // Pre-compute face offsets for optimization
    num_faces_x_ = (ex_ + 1) * ey_ * ez_;
    num_faces_y_ = ex_ * (ey_ + 1) * ez_;
    offset_y_ = num_faces_x_;
    offset_z_ = num_faces_x_ + num_faces_y_;
  }

  PROXY_HOST_DEVICE ModelStruct(const ModelStruct&) = default;
  PROXY_HOST_DEVICE ModelStruct& operator=(const ModelStruct&) = default;
  PROXY_HOST_DEVICE ~ModelStruct() = default;

  /**
   * @brief Convert linear element index to 3D element coordinates
   * @param linearIndex Linear element index [0, ex*ey*ez)
   * @return 3D element index [i, j, k]
   */
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

  /**
   * @brief Get global vertex index from element and local vertex coordinates
   * @param e Element 3D index [i, j, k]
   * @param i Local i-coordinate (0 or 1)
   * @param j Local j-coordinate (0 or 1)
   * @param k Local k-coordinate (0 or 1)
   * @return Global vertex 3D index
   */
  PROXY_HOST_DEVICE
  IndexType globalVertexIndex(IndexType e, int const i, int const j,
                              int const k) const
  {
    return {e[0] + i, e[1] + j, e[2] + k};
  }

  /**
   * @brief Get vertex coordinates from global vertex index
   * @param dofGlobal Global vertex 3D index
   * @param coords Output array for coordinates [x, y, z]
   */
  PROXY_HOST_DEVICE
  void vertexCoords(IndexType dofGlobal, FloatType* const coords) const
  {
    coords[0] = dofGlobal[0] * hx_ + ox_;
    coords[1] = dofGlobal[1] * hy_ + oy_;
    coords[2] = dofGlobal[2] * hz_ + oz_;
  }

  /**
   * @brief Get node coordinate in specified dimension
   *
   * Computes physical coordinate using GLL (Gauss-Lobatto-Legendre) points
   * for spectral element accuracy.
   *
   * @param dofGlobal Global node index (linear)
   * @param dim Dimension (0=x, 1=y, 2=z)
   * @return Physical coordinate value
   */
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

    // offset with local origin
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

  /**
   * @brief Get global node index from element and local coordinates
   * @param e Element linear index
   * @param i Local i-coordinate [0, Order]
   * @param j Local j-coordinate [0, Order]
   * @param k Local k-coordinate [0, Order]
   * @return Global node index (linear)
   */
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

  /**
   * @brief Get P-wave velocity at node
   * @param n Node index
   * @return P-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVpOnNodes(ScalarType n) const
  {
    return 1500;
  }

  /**
   * @brief Get P-wave velocity at element
   * @param e Element index
   * @return P-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVpOnElement(ScalarType e) const
  {
    return 1500;
  }

  /**
   * @brief Get density at node
   * @param n Node index
   * @return Density (kg/m³)
   */
  PROXY_HOST_DEVICE FloatType getModelRhoOnNodes(ScalarType n) const
  {
    return 1;
  }

  /**
   * @brief Get density at element
   * @param e Element index
   * @return Density (kg/m³)
   */
  PROXY_HOST_DEVICE FloatType getModelRhoOnElement(ScalarType e) const
  {
    return 1;
  }

  /**
   * @brief Get S-wave velocity at node
   * @param n Node index
   * @return S-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVsOnNodes(ScalarType n) const
  {
    return 755;
  }

  /**
   * @brief Get S-wave velocity at element
   * @param e Element index
   * @return S-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVsOnElement(ScalarType e) const
  {
    return 755;
  }

  /**
   * @brief Get Thomsen delta parameter at node
   * @param n Node index
   * @return Thomsen delta (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelDeltaOnNodes(ScalarType n) const
  {
    return 0.0;
  }

  /**
   * @brief Get Thomsen delta parameter at element
   * @param e Element index
   * @return Thomsen delta (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelDeltaOnElement(ScalarType e) const
  {
    return 0.0;
  }

  /**
   * @brief Get Thomsen epsilon parameter at node
   * @param n Node index
   * @return Thomsen epsilon (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnNodes(ScalarType n) const
  {
    return 0.0;
  }

  /**
   * @brief Get Thomsen epsilon parameter at element
   * @param e Element index
   * @return Thomsen epsilon (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnElement(ScalarType e) const
  {
    return 0.0;
  }

  /**
   * @brief Get Thomsen gamma parameter at node
   * @param n Node index
   * @return Thomsen gamma (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelGammaOnNodes(ScalarType n) const
  {
    return 0.0;
  }

  /**
   * @brief Get Thomsen gamma parameter at element
   * @param e Element index
   * @return Thomsen gamma (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelGammaOnElement(ScalarType e) const
  {
    return 0.0;
  }

  /**
   * @brief Get tilt angle at node
   * @param n Node index
   * @return Tilt angle theta (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelThetaOnNodes(ScalarType n) const
  {
    return 0;
  }

  /**
   * @brief Get tilt angle at element
   * @param e Element index
   * @return Tilt angle theta (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelThetaOnElement(ScalarType e) const
  {
    return 0;
  }

  /**
   * @brief Get azimuth angle at node
   * @param n Node index
   * @return Azimuth angle phi (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelPhiOnNodes(ScalarType n) const
  {
    return 0.0;
  }

  /**
   * @brief Get azimuth angle at element
   * @param e Element index
   * @return Azimuth angle phi (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelPhiOnElement(ScalarType e) const
  {
    return 0.0;
  }

  /**
   * @brief Initialize elasticity tensors for all elements
   *
   * Computes 6x6 Voigt elasticity tensor from uniform material properties.
   * Only executed if isElastic_ is true.
   */
  void initElasticityTensors()
  {
    if (!isElastic_) return;

    int n_element = ex_ * ey_ * ez_;
    model_C_tensor_element_ = allocateArray3D<array3DReal>(n_element, 6, 6);
    auto& C_tensor = model_C_tensor_element_;

    MAINLOOPHEAD(n_element, i)
    FloatType CTTI[6][6];
    FloatType vp = 1500.0, vs = 755.0, rho = 1.0;
    FloatType delta = 0., epsilon = 0., gamma = 0.0, theta = 0.0, phi = 0.0;
    computeCTensor(vp, vs, rho, delta, epsilon, gamma, theta, phi, CTTI);
    for (int k = 0; k < 6; k++)
      for (int l = 0; l < 6; l++) C_tensor(i, k, l) = CTTI[k][l];
    MAINLOOPEND
  }

  /**
   * @brief Get elasticity tensor for element
   * @param e Element index
   * @param CTTI Output 6x6 Voigt elasticity tensor
   */
  PROXY_HOST_DEVICE
  void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const
  {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) CTTI[i][j] = model_C_tensor_element_(e, i, j);
  }

  /**
   * @brief Get total number of elements
   * @return ex * ey * ez
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfElements() const
  {
    return ex_ * ey_ * ez_;
  }

  /**
   * @brief Get total number of nodes
   * @return (Order*ex+1) * (Order*ey+1) * (Order*ez+1)
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const
  {
    return (Order * ex_ + 1) * (Order * ey_ + 1) * (Order * ez_ + 1);
  }

  /**
   * @brief Get number of nodes per element
   * @return (Order+1)³
   */
  PROXY_HOST_DEVICE int getNumberOfPointsPerElement() const
  {
    return (Order + 1) * (Order + 1) * (Order + 1);
  }

  /**
   * @brief Get polynomial order of elements
   * @return Template parameter Order
   */
  PROXY_HOST_DEVICE int getOrder() const { return Order; }

  /**
   * @brief Compute outward normal vector for element face
   *
   * For Cartesian meshes, normals are axis-aligned: ±[1,0,0], ±[0,1,0],
   * ±[0,0,1]
   *
   * @param e Element index
   * @param local_face Face identifier (kXMinus, kXPlus, etc.)
   * @param v Output normal vector [nx, ny, nz] (normalized)
   */
  PROXY_HOST_DEVICE
  void faceNormal(ScalarType e, CubicFace local_face,
                  FloatType v[3]) const final
  {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;

    int direction = static_cast<int>(local_face) / 2;  // 0=x, 1=y, 2=z
    FloatType sign = (static_cast<int>(local_face) % 2) ? +1.0 : -1.0;

    v[direction] = sign;
  }

  /**
   * @brief Get domain size in specified dimension
   * @param dim Dimension (0=x, 1=y, 2=z)
   * @return Domain size (meters)
   */
  PROXY_HOST_DEVICE FloatType domainSize(int dim) const
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

  /**
   * @brief Compute minimum node spacing in mesh
   *
   * Uses order-dependent coefficients for GLL node spacing.
   *
   * @return Minimum spacing between adjacent GLL nodes (meters)
   */
  PROXY_HOST_DEVICE FloatType getMinSpacing() const
  {
    if constexpr (Order == 1) return min(hx_, min(hy_, hz_));
    if constexpr (Order == 2) return min(hx_, min(hy_, hz_)) / 2;
    if constexpr (Order == 3) return min(hx_, min(hy_, hz_)) * 0.276393;
    if constexpr (Order == 4) return min(hx_, min(hy_, hz_)) * 0.172673;
    if constexpr (Order == 5) return min(hx_, min(hy_, hz_)) * 0.117472;
    return -1;
  }

  /**
   * @brief Get maximum wave speed in mesh
   * @return Uniform P-wave velocity 1500 m/s
   */
  FloatType getMaxSpeed() const { return 1500; }

  /**
   * @brief Check if material properties are stored on nodes
   * @return True if on nodes, false if on elements
   */
  PROXY_HOST_DEVICE bool isModelOnNodes() const { return isModelOnNodes_; }

  /**
   * @brief Check if mesh is for elastic wave propagation
   * @return True if elastic, false if acoustic
   */
  PROXY_HOST_DEVICE bool isElastic() const { return isElastic_; }

  // ============================================================================
  // FACE CONNECTIVITY FUNCTIONS (Computed on-the-fly for Cartesian meshes)
  // ============================================================================

  /**
   * @brief Get global face ID from element and local face (computed on-the-fly)
   *
   * Face numbering scheme for Cartesian meshes:
   * - X-direction faces: [0, (ex+1)*ey*ez)
   * - Y-direction faces: [(ex+1)*ey*ez, (ex+1)*ey*ez + ex*(ey+1)*ez)
   * - Z-direction faces: [remaining...]
   *
   * Uses pre-computed offsets for optimization.
   *
   * @param elem_linear Element linear index
   * @param local_face Local face identifier (0-5)
   * @return Global face ID
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalFace(ScalarType elem_linear, CubicFace local_face) const
  {
    ScalarType ez = elem_linear / (ex_ * ey_);
    ScalarType tmp = elem_linear % (ex_ * ey_);
    ScalarType ey = tmp / ex_;
    ScalarType ex = tmp % ex_;

    switch (local_face)
    {
      case CubicFace::kXMinus:
        return ex + ey * (ex_ + 1) + ez * (ex_ + 1) * ey_;
      case CubicFace::kXPlus:
        return (ex + 1) + ey * (ex_ + 1) + ez * (ex_ + 1) * ey_;
      case CubicFace::kYMinus:
        return offset_y_ + ex + ey * ex_ + ez * ex_ * (ey_ + 1);
      case CubicFace::kYPlus:
        return offset_y_ + ex + (ey + 1) * ex_ + ez * ex_ * (ey_ + 1);
      case CubicFace::kZMinus:
        return offset_z_ + ex + ey * ex_ + ez * ex_ * ey_;
      case CubicFace::kZPlus:
        return offset_z_ + ex + ey * ex_ + (ez + 1) * ex_ * ey_;
      default:
        return -1;
    }
  }

  /**
   * @brief Get global node index from face and local DOF (computed on-the-fly)
   *
   * Uses pre-computed offsets for optimization.
   *
   * @param face_global Global face ID
   * @param local_dof Local node index on face [0, (Order+1)²)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalNodeFromFace(ScalarType face_global, int local_dof) const
  {
    if (face_global < num_faces_x_)
    {
      // X-direction face
      ScalarType face_idx = face_global;
      ScalarType i = face_idx % (ex_ + 1);
      ScalarType j = (face_idx / (ex_ + 1)) % ey_;
      ScalarType k = face_idx / ((ex_ + 1) * ey_);

      int dj = local_dof % (Order + 1);
      int dk = local_dof / (Order + 1);

      ScalarType ni = i * Order;
      ScalarType nj = j * Order + dj;
      ScalarType nk = k * Order + dk;

      return ni + nj * nx_ + nk * nx_ * ny_;
    }
    else if (face_global < offset_z_)
    {
      // Y-direction face
      ScalarType face_idx = face_global - offset_y_;
      ScalarType i = face_idx % ex_;
      ScalarType j = (face_idx / ex_) % (ey_ + 1);
      ScalarType k = face_idx / (ex_ * (ey_ + 1));

      int di = local_dof % (Order + 1);
      int dk = local_dof / (Order + 1);

      ScalarType ni = i * Order + di;
      ScalarType nj = j * Order;
      ScalarType nk = k * Order + dk;

      return ni + nj * nx_ + nk * nx_ * ny_;
    }
    else
    {
      // Z-direction face
      ScalarType face_idx = face_global - offset_z_;
      ScalarType i = face_idx % ex_;
      ScalarType j = (face_idx / ex_) % ey_;
      ScalarType k = face_idx / (ex_ * ey_);

      int di = local_dof % (Order + 1);
      int dj = local_dof / (Order + 1);

      ScalarType ni = i * Order + di;
      ScalarType nj = j * Order + dj;
      ScalarType nk = k * Order;

      return ni + nj * nx_ + nk * nx_ * ny_;
    }
  }

  /**
   * @brief Check if face is on domain boundary (computed on-the-fly)
   *
   * Boundary detection based on Cartesian grid structure:
   * - X-faces: boundary if i=0 or i=ex
   * - Y-faces: boundary if j=0 or j=ey
   * - Z-faces: boundary if k=0 or k=ez
   *
   * Uses pre-computed offsets for optimization.
   *
   * @param face_global Global face ID
   * @return True if boundary face
   */
  PROXY_HOST_DEVICE
  bool isBoundaryFace(ScalarType face_global) const
  {
    if (face_global < num_faces_x_)
    {
      ScalarType i = face_global % (ex_ + 1);
      return (i == 0 || i == ex_);
    }
    else if (face_global < offset_z_)
    {
      ScalarType face_idx = face_global - offset_y_;
      ScalarType j = (face_idx / ex_) % (ey_ + 1);
      return (j == 0 || j == ey_);
    }
    else
    {
      ScalarType face_idx = face_global - offset_z_;
      ScalarType k = face_idx / (ex_ * ey_);
      return (k == 0 || k == ez_);
    }
  }

  /**
   * @brief Get total number of faces in mesh
   *
   * For Cartesian mesh:
   * - X-direction: (ex+1) * ey * ez faces
   * - Y-direction: ex * (ey+1) * ez faces
   * - Z-direction: ex * ey * (ez+1) faces
   *
   * @return Total number of faces
   */
  PROXY_HOST_DEVICE
  ScalarType getNumberOfFaces() const
  {
    return offset_z_ + ex_ * ey_ * (ez_ + 1);
  }

  /**
   * @brief Build face connectivity (no-op for Cartesian meshes)
   *
   * Cartesian meshes compute face connectivity on-the-fly using arithmetic
   * formulas, so no explicit table construction is needed.
   */
  void buildFaceConnectivity() override
  {
    // Nothing to do - faces computed on-the-fly for Cartesian meshes
  }

  /**
   * @brief Check if node is on free surface
   * @param n Node index
   * @return True if free surface node
   */
  PROXY_HOST_DEVICE
  bool isFreeSurface(ScalarType n) const override
  {
    return freeSurfaceTag_[n] == 1;
  }

  /**
   * @brief Enable or disable free surface condition
   * @param enable True to enable free surface, false to disable
   */
  void setFreeSurfaceEnabled(bool enable) override
  {
    freeSurfaceEnabled_ = enable;
  }

  /**
   * @brief Initialize free surface node tags based on mesh geometry
   *
   * Marks nodes as free surface if they are located at the top of the domain
   * (z = oz_ + lz_) within a small tolerance.
   * Considers the freeSurfaceEnabled_ flag to determine if marking is applied.
   *
   * Optimized with parallel loop for GPU execution.
   */
  void initFreeSurface() override
  {
    // Allocate
    freeSurfaceTag_ =
        allocateVector<VECTOR_INT_VIEW>(getNumberOfNodes(), "freeSurfaceTag");

    FloatType tol = getMinSpacing() * 1e-4;
    FloatType z_max = oz_ + lz_;
    bool enabled = freeSurfaceEnabled_;

    // Capture for kernel
    auto tag = freeSurfaceTag_;
    auto mesh_copy = *this;

    LOOPHEAD(getNumberOfNodes(), n)
    {
      FloatType z = mesh_copy.nodeCoord(n, 2);
      bool is_top = (fabs(z - z_max) < tol);
      tag[n] = (is_top && enabled) ? 1 : 0;
    }
    LOOPEND
  }

 private:
  ScalarType ex_, ey_, ez_;
  ScalarType nx_, ny_, nz_;
  FloatType lx_, ly_, lz_;
  FloatType hx_, hy_, hz_;
  FloatType ox_, oy_, oz_;
  bool isModelOnNodes_;
  bool isElastic_;

  // Pre-computed face offsets for optimization
  ScalarType num_faces_x_;
  ScalarType num_faces_y_;
  ScalarType offset_y_;
  ScalarType offset_z_;

  array3DReal model_C_tensor_element_;
  VECTOR_INT_VIEW freeSurfaceTag_;
  bool freeSurfaceEnabled_;
};

}  // namespace model
#endif  // SRC_MODEL_MESH_IMPL_MODEL_STRUCT_MODEL_STRUCT_H_