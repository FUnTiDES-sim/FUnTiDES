#ifndef SRC_MODEL_MODELAPI_INCLUDE_MODEL_UNSTRUCT_H_
#define SRC_MODEL_MODELAPI_INCLUDE_MODEL_UNSTRUCT_H_

#include <elasticity_utils.h>
#include <model.h>

#include <algorithm>
#include <array>
#include <map>
#include <vector>

namespace model
{

// Temporary structure for construction (CPU only)
template <typename ScalarType>
struct FaceDataTemp
{
  std::vector<ScalarType> dofs;
  ScalarType elem_owner = -1;
  ScalarType elem_neighbor = -1;
  int local_face_owner = -1;
  int local_face_neighbor = -1;

  bool isBoundary() const { return elem_neighbor == -1; }
};

// GPU-compatible face connectivity using project macros
template <typename ScalarType>
struct FaceConnectivity
{
  ARRAY_INT_VIEW elem_to_faces_;
  ARRAY_INT_VIEW face_dofs_;
  VECTOR_INT_VIEW face_elem_owner_;
  VECTOR_INT_VIEW face_elem_neighbor_;
  VECTOR_INT_VIEW face_local_owner_;
  VECTOR_INT_VIEW face_local_neighbor_;

  int n_faces_ = 0;
  int ndofs_per_face_ = 0;

  PROXY_HOST_DEVICE
  ScalarType numFaces() const { return n_faces_; }

  /**
   * @brief Map element and local face to global face index
   * @param elem Element index
   * @param local_face Local face number (0-5: 0=x-, 1=x+, 2=y-, 3=y+, 4=z-,
   * 5=z+)
   * @return Global face index
   */
  PROXY_HOST_DEVICE
  ScalarType globalFace(ScalarType elem, int local_face) const
  {
    return elem_to_faces_(elem, local_face);
  }

  /**
   * @brief Map global face and local DOF to global node index
   * @param face_id Global face index
   * @param local_dof Local DOF on face (0 to (order+1)²-1, lexicographic
   * ordering)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType globalNode(ScalarType face_id, int local_dof) const
  {
    return face_dofs_(face_id, local_dof);
  }

  /**
   * @brief Check if a face is on the domain boundary
   * @param face_id Global face index
   * @return True if boundary face (no neighbor element), false if internal
   */
  PROXY_HOST_DEVICE
  bool isBoundary(ScalarType face_id) const
  {
    return face_elem_neighbor_(face_id) == -1;
  }

  /**
   * @brief Get the owner element of a face
   * @param face_id Global face index
   * @return Element index of the owner (first element encountered during
   * construction)
   */
  PROXY_HOST_DEVICE
  ScalarType elemOwner(ScalarType face_id) const
  {
    return face_elem_owner_(face_id);
  }

  /**
   * @brief Get the neighbor element of a face
   * @param face_id Global face index
   * @return Element index of the neighbor, or -1 if boundary face
   */
  PROXY_HOST_DEVICE
  ScalarType elemNeighbor(ScalarType face_id) const
  {
    return face_elem_neighbor_(face_id);
  }

  /**
   * @brief Get the local face number on the owner element
   * @param face_id Global face index
   * @return Local face number (0-5) on the owner element
   */
  PROXY_HOST_DEVICE
  int localFaceOwner(ScalarType face_id) const
  {
    return face_local_owner_(face_id);
  }

  /**
   * @brief Get the local face number on the neighbor element
   * @param face_id Global face index
   * @return Local face number (0-5) on the neighbor element, or -1 if boundary
   */
  PROXY_HOST_DEVICE
  int localFaceNeighbor(ScalarType face_id) const
  {
    return face_local_neighbor_(face_id);
  }
};

template <typename FloatType, typename ScalarType>
struct ModelUnstructData : public ModelDataBase<FloatType, ScalarType>
{
  PROXY_HOST_DEVICE ModelUnstructData() = default;
  PROXY_HOST_DEVICE ~ModelUnstructData() = default;
  PROXY_HOST_DEVICE ModelUnstructData(const ModelUnstructData&) = default;
  PROXY_HOST_DEVICE ModelUnstructData& operator=(const ModelUnstructData&) =
      default;

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
      ARRAY3D_REAL_VIEW model_C_tensor_element, VECTOR_REAL_VIEW boundaries_t)
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
        boundaries_t_(boundaries_t)
  {
  }

  ScalarType order_;
  ScalarType n_element_;
  ScalarType n_node_;
  FloatType lx_, ly_, lz_;
  FloatType ox_, oy_, oz_;
  bool isModelOnNodes_;
  bool isElastic_;

  // Coordinates and index map views
  ARRAY_INT_VIEW global_node_index_;
  VECTOR_REAL_VIEW nodes_coords_x_;
  VECTOR_REAL_VIEW nodes_coords_y_;
  VECTOR_REAL_VIEW nodes_coords_z_;

  // Models view
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
};

/**
 * @brief Unstructured 3D mesh implementation
 */
template <typename FloatType, typename ScalarType>
class ModelUnstruct final : public ModelApi<FloatType, ScalarType>
{
 public:
  /// Define IndexType as an integer for unstructured indexing
  using IndexType = int;

  /**
   * @brief Default constructor.
   */
  PROXY_HOST_DEVICE ModelUnstruct() = default;

  /**
   * @brief Constructor from ModelData.
   * @param data ModelData structure containing all the mesh data
   */
  PROXY_HOST_DEVICE ModelUnstruct(
      const ModelUnstructData<FloatType, ScalarType>& data)
      : order_(data.order_),
        n_element_(data.n_element_),
        n_node_(data.n_node_),
        lx_(data.lx_),
        ly_(data.ly_),
        lz_(data.lz_),
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
        n_points_per_element_((order_ + 1) * (order_ + 1) * (order_ + 1))
  {
  }

  /**
   * @brief Assignment operator.
   */
  PROXY_HOST_DEVICE ModelUnstruct& operator=(const ModelUnstruct&) = default;

  /**
   * @brief Destructor.
   */
  PROXY_HOST_DEVICE ~ModelUnstruct() = default;

  /**
   * @brief pass through function to go from linear index space to IndexType.
   * @param linearIndex The linear index of the element
   * @return the linear index of the element
   */
  PROXY_HOST_DEVICE
  IndexType elementIndex(const int linearIndex) const { return linearIndex; }

  /**
   * @brief Get the global vertex index the element index and local indices.
   * @param e Element index
   * @param i Local i-index of vertex in the element
   * @param j Local j-index of vertex in the element
   * @param k Local k-index of vertex in the element
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
   * @brief Get the vertex coordinates of a global node.
   * @param dofGlobal Global node index
   * @param[out] coords Output array (size 3) holding the coordinates
   */
  PROXY_HOST_DEVICE
  void vertexCoords(IndexType dofGlobal, FloatType* const coords) const
  {
    coords[0] = nodes_coords_x_[dofGlobal];
    coords[1] = nodes_coords_y_[dofGlobal];
    coords[2] = nodes_coords_z_[dofGlobal];
  }

  /**
   * @brief Get the coordinate of a global node in the given dimension.
   * @param dofGlobal Global node index
   * @param dim Dimension index (0 = x, 1 = y, 2 = z)
   * @return Coordinate value in the specified dimension
   */
  PROXY_HOST_DEVICE
  FloatType nodeCoord(ScalarType dofGlobal, int dim) const final
  {
    switch (dim)
    {
      case 0:
        return nodes_coords_x_[dofGlobal];
      case 1:
        return nodes_coords_y_[dofGlobal];
      case 2:
        return nodes_coords_z_[dofGlobal];
      default:
        return FloatType(-1);
    }
  }

  /**
   * @brief Get the global node index for a local element-node triplet.
   * @param e Element index
   * @param i Local i-index in the element
   * @param j Local j-index in the element
   * @param k Local k-index in the element
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType globalNodeIndex(ScalarType e, int i, int j, int k) const final
  {
    const auto localDofIndex =
        i + j * (order_ + 1) + k * (order_ + 1) * (order_ + 1);
    return global_node_index_(e, localDofIndex);
  }

  PROXY_HOST_DEVICE
  FloatType getModelVpOnNodes(ScalarType n) const final
  {
    return model_vp_node_[n];
  }

  PROXY_HOST_DEVICE
  FloatType getModelVpOnElement(ScalarType e) const final
  {
    return model_vp_element_[e];
  }

  PROXY_HOST_DEVICE
  FloatType getModelRhoOnNodes(ScalarType n) const final
  {
    return model_rho_node_[n];
  }

  PROXY_HOST_DEVICE
  FloatType getModelRhoOnElement(ScalarType e) const final
  {
    return model_rho_element_[e];
  }

  PROXY_HOST_DEVICE
  FloatType getModelVsOnNodes(ScalarType n) const final
  {
    return model_vs_node_[n];
  }

  PROXY_HOST_DEVICE
  FloatType getModelVsOnElement(ScalarType e) const final
  {
    return model_vs_element_[e];
  }

  PROXY_HOST_DEVICE
  FloatType getModelDeltaOnNodes(ScalarType n) const final
  {
    return model_delta_node_[n];
  }

  PROXY_HOST_DEVICE
  FloatType getModelDeltaOnElement(ScalarType e) const final
  {
    return model_delta_element_[e];
  }

  PROXY_HOST_DEVICE
  FloatType getModelEpsilonOnNodes(ScalarType n) const final
  {
    return model_epsilon_node_[n];
  }

  PROXY_HOST_DEVICE
  FloatType getModelEpsilonOnElement(ScalarType e) const final
  {
    return model_epsilon_element_[e];
  }

  PROXY_HOST_DEVICE
  FloatType getModelGammaOnNodes(ScalarType n) const final
  {
    return model_gamma_node_[n];
  }

  PROXY_HOST_DEVICE
  FloatType getModelGammaOnElement(ScalarType e) const final
  {
    return model_gamma_element_[e];
  }

  PROXY_HOST_DEVICE
  ScalarType getModelPhiOnNodes(ScalarType n) const final
  {
    return model_phi_node_[n];
  }

  PROXY_HOST_DEVICE
  ScalarType getModelPhiOnElement(ScalarType e) const final
  {
    return model_phi_element_[e];
  }

  PROXY_HOST_DEVICE
  ScalarType getModelThetaOnNodes(ScalarType n) const final
  {
    return model_theta_node_[n];
  }

  PROXY_HOST_DEVICE
  ScalarType getModelThetaOnElement(ScalarType e) const final
  {
    return model_theta_element_[e];
  }

  void initElasticityTensors()
  {
    if (!isElastic_) return;

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

  PROXY_HOST_DEVICE
  void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const final
  {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) CTTI[i][j] = model_C_tensor_element_(e, i, j);
  }

  PROXY_HOST_DEVICE
  bool isModelOnNodes() const final { return isModelOnNodes_; }

  PROXY_HOST_DEVICE
  bool isElastic() const final { return isElastic_; }

  PROXY_HOST_DEVICE
  ScalarType getNumberOfElements() const final { return n_element_; }

  PROXY_HOST_DEVICE
  ScalarType getNumberOfNodes() const final { return n_node_; }

  PROXY_HOST_DEVICE
  int getNumberOfPointsPerElement() const final
  {
    return n_points_per_element_;
  }

  PROXY_HOST_DEVICE
  int getOrder() const final { return static_cast<int>(order_); }

  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(ScalarType n) const final
  {
    return static_cast<BoundaryFlag>(boundaries_t_[n]);
  }

  PROXY_HOST_DEVICE
  void faceNormal(ScalarType e, int local_face, FloatType v[3]) const
  {
    ScalarType n0, n1, n2;
    const int o = order_;

    switch (local_face)
    {
      case 0:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, 0, o, 0);
        n2 = globalNodeIndex(e, 0, 0, o);
        break;
      case 1:
        n0 = globalNodeIndex(e, o, 0, 0);
        n1 = globalNodeIndex(e, o, 0, o);
        n2 = globalNodeIndex(e, o, o, 0);
        break;
      case 2:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, 0, 0, o);
        n2 = globalNodeIndex(e, o, 0, 0);
        break;
      case 3:
        n0 = globalNodeIndex(e, 0, o, 0);
        n1 = globalNodeIndex(e, o, o, 0);
        n2 = globalNodeIndex(e, 0, o, o);
        break;
      case 4:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, o, 0, 0);
        n2 = globalNodeIndex(e, 0, o, 0);
        break;
      case 5:
        n0 = globalNodeIndex(e, 0, 0, o);
        n1 = globalNodeIndex(e, 0, o, o);
        n2 = globalNodeIndex(e, o, 0, o);
        break;
      default:
        v[0] = v[1] = v[2] = 0.0;
        return;
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

  PROXY_HOST_DEVICE
  FloatType domainSize(int dim) const final
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

  PROXY_HOST_DEVICE
  FloatType getMinSpacing() const final
  {
    FloatType minSpacing = std::numeric_limits<FloatType>::max();
    constexpr ScalarType e = 0;

    for (int k = 0; k <= order_; ++k)
    {
      for (int j = 0; j <= order_; ++j)
      {
        for (int i = 0; i < order_; ++i)
        {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i + 1, j, k);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          FloatType spacing = sqrt(dx * dx + dy * dy + dz * dz);
          minSpacing = fmin(minSpacing, spacing);
        }
      }
    }

    for (int k = 0; k <= order_; ++k)
    {
      for (int i = 0; i <= order_; ++i)
      {
        for (int j = 0; j < order_; ++j)
        {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i, j + 1, k);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          FloatType spacing = sqrt(dx * dx + dy * dy + dz * dz);
          minSpacing = fmin(minSpacing, spacing);
        }
      }
    }

    for (int j = 0; j <= order_; ++j)
    {
      for (int i = 0; i <= order_; ++i)
      {
        for (int k = 0; k < order_; ++k)
        {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i, j, k + 1);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          FloatType spacing = sqrt(dx * dx + dy * dy + dz * dz);
          minSpacing = fmin(minSpacing, spacing);
        }
      }
    }

    return minSpacing;
  }

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
   * @brief Build face connectivity map for the unstructured mesh
   * Must be called after construction if face information is needed
   * Constructs two maps:
   * - Map 1: (element, local_face) -> global_face
   * - Map 2: (global_face, local_dof) -> global_node
   */
  void buildFaceConnectivity() override
  {
    // Step 1: Build temporary structures on CPU
    std::vector<FaceDataTemp<ScalarType>> faces_temp;
    std::vector<std::vector<ScalarType>> elem_to_faces_temp(n_element_);

    using FaceKey = std::array<ScalarType, 4>;
    std::map<FaceKey, ScalarType> face_map;

    // Construction loop
    for (ScalarType elem = 0; elem < n_element_; ++elem)
    {
      elem_to_faces_temp[elem].resize(6);

      for (int local_face = 0; local_face < 6; ++local_face)
      {
        auto corners = extractFaceCorners(elem, local_face);
        auto face_key = makeFaceKey(corners);

        auto it = face_map.find(face_key);
        if (it == face_map.end())
        {
          // New face
          ScalarType face_id = faces_temp.size();
          face_map[face_key] = face_id;

          FaceDataTemp<ScalarType> face;
          face.dofs = extractFaceDofs(elem, local_face);
          face.elem_owner = elem;
          face.local_face_owner = local_face;
          faces_temp.push_back(face);

          elem_to_faces_temp[elem][local_face] = face_id;
        }
        else
        {
          // Face already seen (internal face)
          ScalarType face_id = it->second;
          faces_temp[face_id].elem_neighbor = elem;
          faces_temp[face_id].local_face_neighbor = local_face;

          elem_to_faces_temp[elem][local_face] = face_id;
        }
      }
    }

    // Step 2: Allocate and fill Views
    int n_faces = faces_temp.size();
    int ndofs_per_face = (order_ + 1) * (order_ + 1);

    face_connectivity_.n_faces_ = n_faces;
    face_connectivity_.ndofs_per_face_ = ndofs_per_face;

    // Allocate Views
    face_connectivity_.elem_to_faces_ =
        allocateArray2D<ARRAY_INT_VIEW>(n_element_, 6);
    face_connectivity_.face_dofs_ =
        allocateArray2D<ARRAY_INT_VIEW>(n_faces, ndofs_per_face);
    face_connectivity_.face_elem_owner_ =
        allocateVector<VECTOR_INT_VIEW>(n_faces);
    face_connectivity_.face_elem_neighbor_ =
        allocateVector<VECTOR_INT_VIEW>(n_faces);
    face_connectivity_.face_local_owner_ =
        allocateVector<VECTOR_INT_VIEW>(n_faces);
    face_connectivity_.face_local_neighbor_ =
        allocateVector<VECTOR_INT_VIEW>(n_faces);

    // Fill Views directly
    for (ScalarType elem = 0; elem < n_element_; ++elem)
    {
      for (int lf = 0; lf < 6; ++lf)
      {
        face_connectivity_.elem_to_faces_(elem, lf) =
            elem_to_faces_temp[elem][lf];
      }
    }

    for (int face_id = 0; face_id < n_faces; ++face_id)
    {
      const auto& face = faces_temp[face_id];

      face_connectivity_.face_elem_owner_(face_id) = face.elem_owner;
      face_connectivity_.face_elem_neighbor_(face_id) = face.elem_neighbor;
      face_connectivity_.face_local_owner_(face_id) = face.local_face_owner;
      face_connectivity_.face_local_neighbor_(face_id) =
          face.local_face_neighbor;

      for (int dof = 0; dof < ndofs_per_face; ++dof)
      {
        face_connectivity_.face_dofs_(face_id, dof) = face.dofs[dof];
      }
    }
  }

  /**
   * @brief Get global face number from element and local face
   * @param elem Element index
   * @param local_face Local face number (0-5)
   * @return Global face index
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalFace(ScalarType elem, int local_face) const
  {
    return face_connectivity_.globalFace(elem, local_face);
  }

  /**
   * @brief Get global node number from face and local DOF
   * @param face_global Global face index
   * @param local_dof Local DOF index on the face (0 to (order+1)² - 1)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalNodeFromFace(ScalarType face_global, int local_dof) const
  {
    return face_connectivity_.globalNode(face_global, local_dof);
  }

  /**
   * @brief Check if a face is on the boundary
   * @param face_global Global face index
   * @return True if boundary face, false otherwise
   */
  PROXY_HOST_DEVICE
  bool isBoundaryFace(ScalarType face_global) const
  {
    return face_connectivity_.isBoundary(face_global);
  }

  /**
   * @brief Get total number of faces
   * @return Number of faces in the mesh
   */
  PROXY_HOST_DEVICE
  ScalarType getNumberOfFaces() const { return face_connectivity_.numFaces(); }

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

  FaceConnectivity<ScalarType> face_connectivity_;

  /**
   * @brief Extract the 4 corner vertices for a face
   * Face numbering convention:
   * 0: x- (left)   1: x+ (right)
   * 2: y- (front)  3: y+ (back)
   * 4: z- (bottom) 5: z+ (top)
   */
  std::array<ScalarType, 4> extractFaceCorners(ScalarType elem,
                                               int local_face) const
  {
    switch (local_face)
    {
      case 0:  // x-
        return {
            globalVertexIndex(elem, 0, 0, 0), globalVertexIndex(elem, 0, 1, 0),
            globalVertexIndex(elem, 0, 1, 1), globalVertexIndex(elem, 0, 0, 1)};
      case 1:  // x+
        return {
            globalVertexIndex(elem, 1, 0, 0), globalVertexIndex(elem, 1, 1, 0),
            globalVertexIndex(elem, 1, 1, 1), globalVertexIndex(elem, 1, 0, 1)};
      case 2:  // y-
        return {
            globalVertexIndex(elem, 0, 0, 0), globalVertexIndex(elem, 1, 0, 0),
            globalVertexIndex(elem, 1, 0, 1), globalVertexIndex(elem, 0, 0, 1)};
      case 3:  // y+
        return {
            globalVertexIndex(elem, 0, 1, 0), globalVertexIndex(elem, 1, 1, 0),
            globalVertexIndex(elem, 1, 1, 1), globalVertexIndex(elem, 0, 1, 1)};
      case 4:  // z-
        return {
            globalVertexIndex(elem, 0, 0, 0), globalVertexIndex(elem, 1, 0, 0),
            globalVertexIndex(elem, 1, 1, 0), globalVertexIndex(elem, 0, 1, 0)};
      case 5:  // z+
        return {
            globalVertexIndex(elem, 0, 0, 1), globalVertexIndex(elem, 1, 0, 1),
            globalVertexIndex(elem, 1, 1, 1), globalVertexIndex(elem, 0, 1, 1)};
    }
    return {};
  }

  /**
   * @brief Extract ALL DOFs on a face in lexicographic order
   */
  std::vector<ScalarType> extractFaceDofs(ScalarType elem, int local_face) const
  {
    const int o = order_;
    const int ndofs = (o + 1) * (o + 1);
    std::vector<ScalarType> dofs;
    dofs.reserve(ndofs);

    switch (local_face)
    {
      case 0:  // x-
        for (int k = 0; k <= o; ++k)
          for (int j = 0; j <= o; ++j)
            dofs.push_back(globalNodeIndex(elem, 0, j, k));
        break;
      case 1:  // x+
        for (int k = 0; k <= o; ++k)
          for (int j = 0; j <= o; ++j)
            dofs.push_back(globalNodeIndex(elem, o, j, k));
        break;
      case 2:  // y-
        for (int k = 0; k <= o; ++k)
          for (int i = 0; i <= o; ++i)
            dofs.push_back(globalNodeIndex(elem, i, 0, k));
        break;
      case 3:  // y+
        for (int k = 0; k <= o; ++k)
          for (int i = 0; i <= o; ++i)
            dofs.push_back(globalNodeIndex(elem, i, o, k));
        break;
      case 4:  // z-
        for (int j = 0; j <= o; ++j)
          for (int i = 0; i <= o; ++i)
            dofs.push_back(globalNodeIndex(elem, i, j, 0));
        break;
      case 5:  // z+
        for (int j = 0; j <= o; ++j)
          for (int i = 0; i <= o; ++i)
            dofs.push_back(globalNodeIndex(elem, i, j, o));
        break;
    }
    return dofs;
  }

  /**
   * @brief Create canonical face key (sorted vertices for uniqueness)
   */
  std::array<ScalarType, 4> makeFaceKey(std::array<ScalarType, 4> corners) const
  {
    std::sort(corners.begin(), corners.end());
    return corners;
  }
};

}  // namespace model

#endif  // SRC_MODEL_MODELAPI_INCLUDE_MODEL_UNSTRUCT_H_