#ifndef SRC_MODEL_MESH_IMPL_MODEL_UNSTRUCT_MODEL_UNSTRUCT_H_
#define SRC_MODEL_MESH_IMPL_MODEL_UNSTRUCT_MODEL_UNSTRUCT_H_

#include <elasticity_utils.h>
#include <model.h>

#include <algorithm>
#include <array>
#include <map>

#include "face_connectivity_unstruct.h"

namespace model
{

/**
 * @brief Data structure for unstructured mesh initialization
 */
template <typename FloatType, typename ScalarType>
struct ModelUnstructData : public ModelDataBase<FloatType, ScalarType>
{
  PROXY_HOST_DEVICE ModelUnstructData() = default;
  PROXY_HOST_DEVICE ~ModelUnstructData() = default;
  PROXY_HOST_DEVICE ModelUnstructData(const ModelUnstructData&) = default;
  PROXY_HOST_DEVICE ModelUnstructData& operator=(const ModelUnstructData&) =
      default;

  /**
   * @brief Full constructor with all mesh data
   */
  PROXY_HOST_DEVICE
  ModelUnstructData(
      ScalarType order, ScalarType n_element, ScalarType n_node, FloatType lx,
      FloatType ly, FloatType lz, bool isModelOnNodes, bool isElastic,
      ARRAY_INT_VIEW global_node_index, VECTOR_REAL_VIEW nodes_coords_x,
      VECTOR_REAL_VIEW nodes_coords_y, VECTOR_REAL_VIEW nodes_coords_z,
      VECTOR_REAL_VIEW model_vp_node, VECTOR_REAL_VIEW model_vp_element,
      VECTOR_REAL_VIEW model_rho_node, VECTOR_REAL_VIEW model_rho_element,
      VECTOR_REAL_VIEW model_vs_node, VECTOR_REAL_VIEW model_vs_element,
      VECTOR_REAL_VIEW model_delta_node, VECTOR_REAL_VIEW model_delta_element,
      VECTOR_REAL_VIEW model_epsilon_node,
      VECTOR_REAL_VIEW model_epsilon_element, VECTOR_REAL_VIEW model_gamma_node,
      VECTOR_REAL_VIEW model_gamma_element, VECTOR_REAL_VIEW model_theta_node,
      VECTOR_REAL_VIEW model_theta_element, VECTOR_REAL_VIEW model_phi_node,
      VECTOR_REAL_VIEW model_phi_element,
      ARRAY3D_REAL_VIEW model_C_tensor_element, VECTOR_REAL_VIEW boundaries_t,
      FaceConnectivityUnstructData<FloatType, ScalarType> face_connectivity =
          {})
      : order_(order),
        n_element_(n_element),
        n_node_(n_node),
        lx_(lx),
        ly_(ly),
        lz_(lz),
        isModelOnNodes_(isModelOnNodes),
        isElastic_(isElastic),
        global_node_index_(global_node_index),
        nodes_coords_x_(nodes_coords_x),
        nodes_coords_y_(nodes_coords_y),
        nodes_coords_z_(nodes_coords_z),
        model_vp_node_(model_vp_node),
        model_vp_element_(model_vp_element),
        model_rho_node_(model_rho_node),
        model_rho_element_(model_rho_element),
        model_vs_node_(model_vs_node),
        model_vs_element_(model_vs_element),
        model_delta_node_(model_delta_node),
        model_delta_element_(model_delta_element),
        model_epsilon_node_(model_epsilon_node),
        model_epsilon_element_(model_epsilon_element),
        model_gamma_node_(model_gamma_node),
        model_gamma_element_(model_gamma_element),
        model_theta_node_(model_theta_node),
        model_theta_element_(model_theta_element),
        model_phi_node_(model_phi_node),
        model_phi_element_(model_phi_element),
        model_C_tensor_element_(model_C_tensor_element),
        boundaries_t_(boundaries_t),
        face_connectivity_(face_connectivity)
  {
  }

  FloatType origin_x_{0}, origin_y_{0}, origin_z_{0};
  FloatType ox_, oy_, oz_;  // Local origin
  ScalarType order_;
  ScalarType n_element_;
  ScalarType n_node_;
  FloatType lx_, ly_, lz_;  // Local dimensions

  bool isModelOnNodes_;
  bool isElastic_;

  ARRAY_INT_VIEW global_node_index_;
  VECTOR_REAL_VIEW nodes_coords_x_;
  VECTOR_REAL_VIEW nodes_coords_y_;
  VECTOR_REAL_VIEW nodes_coords_z_;

  VECTOR_REAL_VIEW model_vp_node_;
  VECTOR_REAL_VIEW model_vp_element_;
  VECTOR_REAL_VIEW model_rho_node_;
  VECTOR_REAL_VIEW model_rho_element_;
  VECTOR_REAL_VIEW model_vs_node_;
  VECTOR_REAL_VIEW model_vs_element_;
  VECTOR_REAL_VIEW model_delta_node_;
  VECTOR_REAL_VIEW model_delta_element_;
  VECTOR_REAL_VIEW model_epsilon_node_;
  VECTOR_REAL_VIEW model_epsilon_element_;
  VECTOR_REAL_VIEW model_gamma_node_;
  VECTOR_REAL_VIEW model_gamma_element_;
  VECTOR_REAL_VIEW model_theta_node_;
  VECTOR_REAL_VIEW model_theta_element_;
  VECTOR_REAL_VIEW model_phi_node_;
  VECTOR_REAL_VIEW model_phi_element_;
  ARRAY3D_REAL_VIEW model_C_tensor_element_;
  VECTOR_REAL_VIEW boundaries_t_;
  FaceConnectivityUnstructData<FloatType, ScalarType> face_connectivity_;
};

/**
 * @brief Unstructured 3D hexahedral mesh implementation
 */
template <typename FloatType, typename ScalarType>
class ModelUnstruct final : public ModelApi<FloatType, ScalarType>
{
 public:
  using IndexType = int;

  PROXY_HOST_DEVICE ModelUnstruct() = default;

  /**
   * @brief Construct from data structure
   */
  PROXY_HOST_DEVICE ModelUnstruct(
      const ModelUnstructData<FloatType, ScalarType>& data)
      : order_(data.order_),
        n_element_(data.n_element_),
        n_node_(data.n_node_),
        lx_(data.lx_),
        ly_(data.ly_),
        lz_(data.lz_),
        ox_(data.ox_),
        oy_(data.oy_),
        oz_(data.oz_),
        isModelOnNodes_(data.isModelOnNodes_),
        isElastic_(data.isElastic_),
        global_node_index_(data.global_node_index_),
        nodes_coords_x_(data.nodes_coords_x_),
        nodes_coords_y_(data.nodes_coords_y_),
        nodes_coords_z_(data.nodes_coords_z_),
        model_vp_node_(data.model_vp_node_),
        model_vp_element_(data.model_vp_element_),
        model_rho_node_(data.model_rho_node_),
        model_rho_element_(data.model_rho_element_),
        model_vs_node_(data.model_vs_node_),
        model_vs_element_(data.model_vs_element_),
        model_delta_node_(data.model_delta_node_),
        model_delta_element_(data.model_delta_element_),
        model_epsilon_node_(data.model_epsilon_node_),
        model_epsilon_element_(data.model_epsilon_element_),
        model_gamma_node_(data.model_gamma_node_),
        model_gamma_element_(data.model_gamma_element_),
        model_phi_node_(data.model_phi_node_),
        model_phi_element_(data.model_phi_element_),
        model_theta_node_(data.model_theta_node_),
        model_theta_element_(data.model_theta_element_),
        model_C_tensor_element_(data.model_C_tensor_element_),
        boundaries_t_(data.boundaries_t_),
        face_connectivity_(data.face_connectivity_),
        n_points_per_element_((order_ + 1) * (order_ + 1) * (order_ + 1))
  {
  }

  PROXY_HOST_DEVICE ModelUnstruct& operator=(const ModelUnstruct&) = default;
  PROXY_HOST_DEVICE ~ModelUnstruct() = default;

  /**
   * @brief Convert linear element index to element identifier
   * @param linearIndex Linear element index
   * @return Element index (identity for unstructured)
   */
  PROXY_HOST_DEVICE
  IndexType elementIndex(const int linearIndex) const { return linearIndex; }

  /**
   * @brief Get global vertex index from element and local vertex coordinates
   * @param e Element index
   * @param i Local i-coordinate (0 or 1)
   * @param j Local j-coordinate (0 or 1)
   * @param k Local k-coordinate (0 or 1)
   * @return Global vertex index
   */
  PROXY_HOST_DEVICE
  IndexType globalVertexIndex(IndexType e, int const i, int const j,
                              int const k) const
  {
    int local_i = i * order_;
    int local_j = j * order_;
    int local_k = k * order_;
    const auto localDofIndex = local_i + local_j * (order_ + 1) +
                               local_k * (order_ + 1) * (order_ + 1);
    return global_node_index_(e, localDofIndex);
  }

  /**
   * @brief Get vertex coordinates
   * @param dofGlobal Global vertex index
   * @param coords Output array for coordinates [x, y, z]
   */
  PROXY_HOST_DEVICE
  void vertexCoords(IndexType dofGlobal, FloatType* const coords) const
  {
    coords[0] = nodes_coords_x_[dofGlobal];
    coords[1] = nodes_coords_y_[dofGlobal];
    coords[2] = nodes_coords_z_[dofGlobal];
  }

  /**
   * @brief Get node coordinate in specified dimension
   * @param dofGlobal Global node index
   * @param dim Dimension (0=x, 1=y, 2=z)
   * @return Coordinate value
   */
  PROXY_HOST_DEVICE
  FloatType nodeCoord(ScalarType dofGlobal, int dim) const final
  {
    switch (dim)
    {
      case 0: {
        return nodes_coords_x_[dofGlobal];
      }
      case 1: {
        return nodes_coords_y_[dofGlobal];
      }
      case 2: {
        return nodes_coords_z_[dofGlobal];
      }
      default:
        return FloatType(-1);
    }
  }

  /**
   * @brief Get global node index from element and local coordinates
   * @param e Element index
   * @param i Local i-coordinate [0, order]
   * @param j Local j-coordinate [0, order]
   * @param k Local k-coordinate [0, order]
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType globalNodeIndex(ScalarType e, int i, int j, int k) const final
  {
    const auto localDofIndex =
        i + j * (order_ + 1) + k * (order_ + 1) * (order_ + 1);
    return global_node_index_(e, localDofIndex);
  }

  /**
   * @brief Get P-wave velocity at node
   * @param n Node index
   * @return P-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVpOnNodes(ScalarType n) const final
  {
    return model_vp_node_[n];
  }

  /**
   * @brief Get P-wave velocity at element
   * @param e Element index
   * @return P-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVpOnElement(ScalarType e) const final
  {
    return model_vp_element_[e];
  }

  /**
   * @brief Get density at node
   * @param n Node index
   * @return Density (kg/m³)
   */
  PROXY_HOST_DEVICE FloatType getModelRhoOnNodes(ScalarType n) const final
  {
    return model_rho_node_[n];
  }

  /**
   * @brief Get density at element
   * @param e Element index
   * @return Density (kg/m³)
   */
  PROXY_HOST_DEVICE FloatType getModelRhoOnElement(ScalarType e) const final
  {
    return model_rho_element_[e];
  }

  /**
   * @brief Get S-wave velocity at node
   * @param n Node index
   * @return S-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVsOnNodes(ScalarType n) const final
  {
    return model_vs_node_[n];
  }

  /**
   * @brief Get S-wave velocity at element
   * @param e Element index
   * @return S-wave velocity (m/s)
   */
  PROXY_HOST_DEVICE FloatType getModelVsOnElement(ScalarType e) const final
  {
    return model_vs_element_[e];
  }

  /**
   * @brief Get Thomsen delta parameter at node
   * @param n Node index
   * @return Thomsen delta (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelDeltaOnNodes(ScalarType n) const final
  {
    return model_delta_node_[n];
  }

  /**
   * @brief Get Thomsen delta parameter at element
   * @param e Element index
   * @return Thomsen delta (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelDeltaOnElement(ScalarType e) const final
  {
    return model_delta_element_[e];
  }

  /**
   * @brief Get Thomsen epsilon parameter at node
   * @param n Node index
   * @return Thomsen epsilon (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnNodes(ScalarType n) const final
  {
    return model_epsilon_node_[n];
  }

  /**
   * @brief Get Thomsen epsilon parameter at element
   * @param e Element index
   * @return Thomsen epsilon (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnElement(ScalarType e) const final
  {
    return model_epsilon_element_[e];
  }

  /**
   * @brief Get Thomsen gamma parameter at node
   * @param n Node index
   * @return Thomsen gamma (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelGammaOnNodes(ScalarType n) const final
  {
    return model_gamma_node_[n];
  }

  /**
   * @brief Get Thomsen gamma parameter at element
   * @param e Element index
   * @return Thomsen gamma (dimensionless)
   */
  PROXY_HOST_DEVICE FloatType getModelGammaOnElement(ScalarType e) const final
  {
    return model_gamma_element_[e];
  }

  /**
   * @brief Get azimuth angle at node
   * @param n Node index
   * @return Azimuth angle phi (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelPhiOnNodes(ScalarType n) const final
  {
    return model_phi_node_[n];
  }

  /**
   * @brief Get azimuth angle at element
   * @param e Element index
   * @return Azimuth angle phi (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelPhiOnElement(ScalarType e) const final
  {
    return model_phi_element_[e];
  }

  /**
   * @brief Get tilt angle at node
   * @param n Node index
   * @return Tilt angle theta (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelThetaOnNodes(ScalarType n) const final
  {
    return model_theta_node_[n];
  }

  /**
   * @brief Get tilt angle at element
   * @param e Element index
   * @return Tilt angle theta (radians)
   */
  PROXY_HOST_DEVICE ScalarType getModelThetaOnElement(ScalarType e) const final
  {
    return model_theta_element_[e];
  }

  /**
   * @brief Initialize elasticity tensors for all elements
   *
   * Computes the 6x6 Voigt elasticity tensor from material properties.
   * Only executed if isElastic_ is true.
   */
  void initElasticityTensors(AnisotropyType anisotropy_type) final
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
      model_C_tensor_element_ = allocateArray3D<array3DReal>(n_element_, 6, 6);

      auto& C_tensor = model_C_tensor_element_;
      auto& vp = model_vp_element_;
      auto& vs = model_vs_element_;
      auto& rho = model_rho_element_;
      auto& delta = model_delta_element_;
      auto& epsilon = model_epsilon_element_;
      auto& gamma = model_gamma_element_;
      auto& theta = model_theta_element_;
      auto& phi = model_phi_element_;

      MAINLOOPHEAD(n_element_, i)
      FloatType CTTI[6][6];
      FloatType vp_val = static_cast<FloatType>(vp[i]);
      FloatType vs_val = static_cast<FloatType>(vs[i]);
      FloatType rho_val = static_cast<FloatType>(rho[i]);
      FloatType delta_val = static_cast<FloatType>(delta[i]);
      FloatType epsilon_val = static_cast<FloatType>(epsilon[i]);
      FloatType gamma_val = static_cast<FloatType>(gamma[i]);
      FloatType theta_val = static_cast<FloatType>(theta[i]);
      FloatType phi_val = static_cast<FloatType>(phi[i]);

      computeCTensor(vp_val, vs_val, rho_val, delta_val, epsilon_val, gamma_val,
                     theta_val, phi_val, CTTI);

      for (int k = 0; k < 6; k++)
        for (int l = 0; l < 6; l++) C_tensor(i, k, l) = CTTI[k][l];
      MAINLOOPEND
    }
  }

  /**
   * @brief Get elasticity tensor for element
   * @param e Element index
   * @param CTTI Output 6x6 Voigt elasticity tensor
   */
  PROXY_HOST_DEVICE
  void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const final
  {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) CTTI[i][j] = model_C_tensor_element_(e, i, j);
  }

  /**
   * @brief Check if material properties are stored on nodes
   * @return True if on nodes, false if on elements
   */
  PROXY_HOST_DEVICE bool isModelOnNodes() const final
  {
    return isModelOnNodes_;
  }

  /**
   * @brief Check if mesh is for elastic wave propagation
   * @return True if elastic, false if acoustic
   */
  PROXY_HOST_DEVICE bool isElastic() const final { return isElastic_; }

  /**
   * @brief Get total number of elements
   * @return Number of elements
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfElements() const final
  {
    return n_element_;
  }

  /**
   * @brief Get total number of nodes
   * @return Number of nodes
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const final
  {
    return n_node_;
  }

  /**
   * @brief Get number of nodes per element
   * @return (order+1)³ nodes per element
   */
  PROXY_HOST_DEVICE int getNumberOfPointsPerElement() const final
  {
    return n_points_per_element_;
  }

  /**
   * @brief Get polynomial order of elements
   * @return Element order
   */
  PROXY_HOST_DEVICE int getOrder() const final
  {
    return static_cast<int>(order_);
  }

  /**
   * @brief Compute outward normal vector for element face
   * @param e Element index
   * @param local_face Face identifier (kXMinus, kXPlus, etc.)
   * @param v Output normal vector [nx, ny, nz] (normalized)
   */
  PROXY_HOST_DEVICE
  void faceNormal(ScalarType e, CubicFace local_face,
                  FloatType v[3]) const final
  {
    ScalarType n0, n1, n2;
    const int o = order_;

    switch (local_face)
    {
      case CubicFace::kXMinus:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, 0, o, 0);
        n2 = globalNodeIndex(e, 0, 0, o);
        break;
      case CubicFace::kXPlus:
        n0 = globalNodeIndex(e, o, 0, 0);
        n1 = globalNodeIndex(e, o, 0, o);
        n2 = globalNodeIndex(e, o, o, 0);
        break;
      case CubicFace::kYMinus:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, 0, 0, o);
        n2 = globalNodeIndex(e, o, 0, 0);
        break;
      case CubicFace::kYPlus:
        n0 = globalNodeIndex(e, 0, o, 0);
        n1 = globalNodeIndex(e, o, o, 0);
        n2 = globalNodeIndex(e, 0, o, o);
        break;
      case CubicFace::kZMinus:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, o, 0, 0);
        n2 = globalNodeIndex(e, 0, o, 0);
        break;
      case CubicFace::kZPlus:
        n0 = globalNodeIndex(e, 0, 0, o);
        n1 = globalNodeIndex(e, 0, o, o);
        n2 = globalNodeIndex(e, o, 0, o);
        break;
    }

    FloatType p0[3], p1[3], p2[3];
    for (int d = 0; d < 3; ++d)
    {
      p0[d] = nodeCoord(n0, d);
      p1[d] = nodeCoord(n1, d);
      p2[d] = nodeCoord(n2, d);
    }

    FloatType t1[3], t2[3];
    t1[0] = p1[0] - p0[0];
    t1[1] = p1[1] - p0[1];
    t1[2] = p1[2] - p0[2];
    t2[0] = p2[0] - p0[0];
    t2[1] = p2[1] - p0[1];
    t2[2] = p2[2] - p0[2];

    v[0] = t1[1] * t2[2] - t1[2] * t2[1];
    v[1] = t1[2] * t2[0] - t1[0] * t2[2];
    v[2] = t1[0] * t2[1] - t1[1] * t2[0];

    FloatType norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm > 1e-12)
    {
      v[0] /= norm;
      v[1] /= norm;
      v[2] /= norm;
    }
  }

  /**
   * @brief Get boundary type flag for a node
   * @param n Node index
   * @return BoundaryFlag enum value
   */
  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(ScalarType n) const override
  {
    if (boundaries_t_.extent(0) == 0) return BoundaryFlag::InteriorNode;
    return static_cast<BoundaryFlag>(static_cast<uint8_t>(boundaries_t_[n]));
  }

  /**
   * @brief Initialize boundary flags based on node positions
   *
   * Detects boundary nodes using geometry and marks them:
   * - Nodes at global domain edges are boundary nodes
   * - Top surface (Z+) marked as Surface if free_surface_on_top=true
   * - Other boundaries marked as Damping (absorbing boundary)
   * - Interior nodes (including MPI inter-domain boundaries) marked as
   * InteriorNode
   *
   * This geometric detection is MPI-safe: only nodes at the GLOBAL domain
   * boundaries are marked, not nodes at MPI partition boundaries.
   *
   * @param free_surface_on_top If true, mark top (Z+) as Surface, else as
   * Damping
   */
  void initializeBoundaryFlags(bool free_surface_on_top) override
  {
    if (boundaries_t_.extent(0) > 0) return;

    boundaries_t_ = allocateVector<VECTOR_REAL_VIEW>(n_node_, "boundaries_t");

    FloatType tol = getMinSpacing() * 1e-4;
    FloatType x_min = ox_, x_max = ox_ + lx_;
    FloatType y_min = oy_, y_max = oy_ + ly_;
    FloatType z_min = oz_, z_max = oz_ + lz_;
    bool enabled_fs = free_surface_on_top;

    auto boundaries = boundaries_t_;
    auto mesh_copy = *this;

    LOOPHEAD(n_node_, n)
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
        boundaries[n] = static_cast<FloatType>(BoundaryFlag::InteriorNode);
      else if (at_zmax && enabled_fs)
        boundaries[n] = static_cast<FloatType>(BoundaryFlag::Surface);
      else
        boundaries[n] = static_cast<FloatType>(BoundaryFlag::Damping);
    }
    LOOPEND
  }

  /**
   * @brief Get domain size in specified dimension
   * @param dim Dimension (0=x, 1=y, 2=z)
   * @return Domain size (meters)
   */
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
        return FloatType(-1);
    }
  }

  /**
   * @brief Compute minimum node spacing in mesh
   * @return Minimum spacing between adjacent nodes (meters)
   */
  PROXY_HOST_DEVICE FloatType getMinSpacing() const final
  {
    FloatType minSpacing = std::numeric_limits<FloatType>::max();
    constexpr ScalarType e = 0;

    for (int k = 0; k <= order_; ++k)
      for (int j = 0; j <= order_; ++j)
        for (int i = 0; i < order_; ++i)
        {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i + 1, j, k);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          minSpacing = fmin(minSpacing, sqrt(dx * dx + dy * dy + dz * dz));
        }

    for (int k = 0; k <= order_; ++k)
      for (int i = 0; i <= order_; ++i)
        for (int j = 0; j < order_; ++j)
        {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i, j + 1, k);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          minSpacing = fmin(minSpacing, sqrt(dx * dx + dy * dy + dz * dz));
        }

    for (int j = 0; j <= order_; ++j)
      for (int i = 0; i <= order_; ++i)
        for (int k = 0; k < order_; ++k)
        {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i, j, k + 1);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          minSpacing = fmin(minSpacing, sqrt(dx * dx + dy * dy + dz * dz));
        }

    return minSpacing;
  }

  /**
   * @brief Find maximum wave speed in mesh
   * @return Maximum P-wave velocity (m/s)
   */
  FloatType getMaxSpeed() const final
  {
    FloatType maxSpeedNode = std::numeric_limits<FloatType>::lowest();
    FloatType maxSpeedElem = std::numeric_limits<FloatType>::lowest();

    if (model_vp_node_.extent(0) > 0)
    {
      FIND_MAX_1D(model_vp_node_, n_node_, maxSpeedNode);
    }
    else if (model_vp_element_.extent(0) > 0)
    {
      FIND_MAX_1D(model_vp_element_, n_element_, maxSpeedElem);
    }
    else
    {
      throw std::runtime_error(
          "No model initialized (model unstruct getMaxSpeed).");
    }
    return max(maxSpeedElem, maxSpeedNode);
  }

  // ============================================================================
  // FACE CONNECTIVITY FUNCTIONS
  // ============================================================================

  /**
   * @brief Build face connectivity tables for absorbing boundary conditions
   *
   * Constructs mappings between elements, faces, and nodes:
   * - Identifies unique faces across all elements
   * - Detects shared internal faces between adjacent elements
   * - Flags boundary faces (faces with only one adjacent element)
   * - Builds node lists for each face (for quadrature integration)
   *
   * Must be called before using face-related query functions.
   * Complexity: O(N_elem) with map-based face matching.
   */
  void buildFaceConnectivity() override
  {
    if (face_connectivity_.getNumberOfFaces() > 0) return;

    face_connectivity_.build(*this);
  }
  /**
   * @brief Get global face ID from element and local face
   * @param elem Element index
   * @param local_face Local face identifier (0-5)
   * @return Global face ID
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const
  {
    return face_connectivity_.getGlobalFace(elem, local_face);
  }

  /**
   * @brief Get global node index from face and local DOF
   * @param face_global Global face ID
   * @param local_dof Local node index on face [0, (order+1)²)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalNodeFromFace(ScalarType face_global, int local_dof) const
  {
    return face_connectivity_.getGlobalNodeFromFace(face_global, local_dof);
  }

  /**
   * @brief Check if face is on domain boundary
   * @param face_global Global face ID
   * @return True if boundary face (no neighbor element)
   */
  PROXY_HOST_DEVICE
  bool isBoundaryFace(ScalarType face_global) const
  {
    return face_connectivity_.isBoundaryFace(face_global);
  }

  /**
   * @brief Get total number of faces in mesh
   * @return Number of unique faces
   */
  PROXY_HOST_DEVICE
  ScalarType getNumberOfFaces() const
  {
    return face_connectivity_.getNumberOfFaces();
  }

  /**
   * @brief Check if node is on free surface
   * @param n Node index
   * @return True if free surface node, false otherwise
   */
  PROXY_HOST_DEVICE
  bool isFreeSurface(ScalarType n) const override
  {
    return freeSurfaceTag_[n] == 1;
  }

  /**
   * @brief Enable or disable free surface marking
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
   */

  void initFreeSurface() override
  {
    if (freeSurfaceTag_.extent(0) > 0) return;

    freeSurfaceTag_ =
        allocateVector<VECTOR_INT_VIEW>(getNumberOfNodes(), "freeSurfaceTag");

    FloatType tol = getMinSpacing() * 1e-4;
    FloatType z_max = oz_ + lz_;
    bool enabled = freeSurfaceEnabled_;

    auto tag = freeSurfaceTag_;
    auto mesh_copy = *this;

    LOOPHEAD(getNumberOfNodes(), n)
    {
      FloatType z = mesh_copy.nodeCoord(n, 2);
      tag[n] = (fabs(z - z_max) < tol && enabled) ? 1 : 0;
    }
    LOOPEND
  }

 private:
  ScalarType order_;
  ScalarType n_element_;
  ScalarType n_node_;
  FloatType lx_, ly_, lz_;
  FloatType ox_, oy_, oz_;

  int n_points_per_element_;
  bool isModelOnNodes_;
  bool isElastic_;

  ARRAY_INT_VIEW global_node_index_;
  VECTOR_REAL_VIEW nodes_coords_x_;
  VECTOR_REAL_VIEW nodes_coords_y_;
  VECTOR_REAL_VIEW nodes_coords_z_;

  VECTOR_REAL_VIEW model_vp_node_;
  VECTOR_REAL_VIEW model_vp_element_;
  VECTOR_REAL_VIEW model_rho_node_;
  VECTOR_REAL_VIEW model_rho_element_;
  VECTOR_REAL_VIEW model_vs_node_;
  VECTOR_REAL_VIEW model_vs_element_;
  VECTOR_REAL_VIEW model_delta_node_;
  VECTOR_REAL_VIEW model_delta_element_;
  VECTOR_REAL_VIEW model_epsilon_node_;
  VECTOR_REAL_VIEW model_epsilon_element_;
  VECTOR_REAL_VIEW model_theta_node_;
  VECTOR_REAL_VIEW model_theta_element_;
  VECTOR_REAL_VIEW model_gamma_node_;
  VECTOR_REAL_VIEW model_gamma_element_;
  VECTOR_REAL_VIEW model_phi_node_;
  VECTOR_REAL_VIEW model_phi_element_;
  ARRAY3D_REAL_VIEW model_C_tensor_element_;
  VECTOR_REAL_VIEW boundaries_t_;
  VECTOR_INT_VIEW freeSurfaceTag_;
  bool freeSurfaceEnabled_;

  FaceConnectivityUnstruct<FloatType, ScalarType> face_connectivity_;
};

}  // namespace model

#endif  // SRC_MODEL_MESH_IMPL_MODEL_UNSTRUCT_MODEL_UNSTRUCT_H_