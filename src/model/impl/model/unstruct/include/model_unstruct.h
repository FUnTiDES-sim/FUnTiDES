#ifndef SRC_MODEL_MODELAPI_INCLUDE_MODEL_UNSTRUCT_H_
#define SRC_MODEL_MODELAPI_INCLUDE_MODEL_UNSTRUCT_H_

#include <model.h>

namespace model
{
template <typename FloatType, typename ScalarType>
struct ModelUnstructData : public ModelDataBase<FloatType, ScalarType>
{
  // GPU-compatible special member functions
  PROXY_HOST_DEVICE ModelUnstructData() = default;
  PROXY_HOST_DEVICE ~ModelUnstructData() = default;
  PROXY_HOST_DEVICE ModelUnstructData(const ModelUnstructData&) = default;
  PROXY_HOST_DEVICE ModelUnstructData& operator=(const ModelUnstructData&) =
      default;

  PROXY_HOST_DEVICE
  ModelUnstructData(
      ScalarType order, ScalarType n_element, ScalarType n_node, FloatType lx,
      FloatType ly, FloatType lz, bool isModelOnNodes,
      ARRAY_INT_VIEW global_node_index, VECTOR_REAL_VIEW nodes_coords_x,
      VECTOR_REAL_VIEW nodes_coords_y, VECTOR_REAL_VIEW nodes_coords_z,
      VECTOR_REAL_VIEW model_vp_node, VECTOR_REAL_VIEW model_vp_element,
      VECTOR_REAL_VIEW model_rho_node, VECTOR_REAL_VIEW model_rho_element,
      VECTOR_REAL_VIEW boundaries_t)
      : order_(order),
        n_element_(n_element),
        n_node_(n_node),
        lx_(lx),
        ly_(ly),
        lz_(lz),
        isModelOnNodes_(isModelOnNodes),
        global_node_index_(global_node_index),
        nodes_coords_x_(nodes_coords_x),
        nodes_coords_y_(nodes_coords_y),
        nodes_coords_z_(nodes_coords_z),
        model_vp_node_(model_vp_node),
        model_vp_element_(model_vp_element),
        model_rho_node_(model_rho_node),
        model_rho_element_(model_rho_element),
        boundaries_t_(boundaries_t)
  {
  }

  ScalarType order_;
  ScalarType n_element_;
  ScalarType n_node_;
  FloatType lx_, ly_, lz_;
  bool isModelOnNodes_;

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
  VECTOR_REAL_VIEW boundaries_t_;
};

/**
 * @brief Abstract base class representing a structured 3D mesh.
 */
template <typename FloatType, typename ScalarType>
class ModelUnstruct : public ModelApi<FloatType, ScalarType>
{
 public:
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
        global_node_index_(data.global_node_index_),
        nodes_coords_x_(data.nodes_coords_x_),
        nodes_coords_y_(data.nodes_coords_y_),
        nodes_coords_z_(data.nodes_coords_z_),
        model_vp_node_(data.model_vp_node_),
        model_vp_element_(data.model_vp_element_),
        model_rho_node_(data.model_rho_node_),
        model_rho_element_(data.model_rho_element_),
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
      case 0: {
        return nodes_coords_x_[dofGlobal];  // Fixed: was dofGlobalIndex
      }
      case 1: {
        return nodes_coords_y_[dofGlobal];  // Fixed: was dofGlobalIndex
      }
      case 2: {
        return nodes_coords_z_[dofGlobal];  // Fixed: was dofGlobalIndex
      }
      default:
        return FloatType(-1);  // Cast to proper type
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
    return global_node_index_(e, localDofIndex);  // Fixed: was elementIndex
  }

  /**
   * @brief Get the P-wave velocity value at a global node.
   * @param n Global node index
   * @return Model P-wave velocity value at the node
   */
  PROXY_HOST_DEVICE
  FloatType getModelVpOnNodes(ScalarType n) const final
  {  // Added const, removed virtual
    return model_vp_node_[n];
  }

  /**
   * @brief Get the average P-wave velocity value on a given element.
   * @param e Element index
   * @return Model P-wave velocity value for the element
   */
  PROXY_HOST_DEVICE
  FloatType getModelVpOnElement(ScalarType e) const final
  {                               // Added const, removed virtual, fixed typo
    return model_vp_element_[e];  // Fixed: was "reutrn"
  }

  /**
   * @brief Get the density value at a global node.
   * @param n Global node index
   * @return Model density value at the node
   */
  PROXY_HOST_DEVICE
  FloatType getModelRhoOnNodes(ScalarType n) const final
  {  // Added const, removed virtual
    return model_rho_node_[n];
  }

  /**
   * @brief Get the average density value on a given element.
   * @param e Element index
   * @return Model density value for the element
   */
  PROXY_HOST_DEVICE
  FloatType getModelRhoOnElement(ScalarType e) const final
  {  // Added const, removed virtual
    return model_rho_element_[e];
  }

  PROXY_HOST_DEVICE
  bool isModelOnNodes() const final { return isModelOnNodes_; }

  /**
   * @brief Get the total number of elements in the mesh.
   * @return Total element count
   */
  PROXY_HOST_DEVICE
  ScalarType getNumberOfElements() const final
  {  // Added const, removed virtual
    return n_element_;
  }

  /**
   * @brief Get the total number of global nodes in the mesh.
   * @return Total node count
   */
  PROXY_HOST_DEVICE
  ScalarType getNumberOfNodes() const final
  {  // Added const, removed virtual
    return n_node_;
  }

  /**
   * @brief Get the number of interpolation points per element.
   * @return Number of interpolation points in one element
   */
  PROXY_HOST_DEVICE
  int getNumberOfPointsPerElement() const final
  {
    return n_points_per_element_;
  }

  /**
   * @brief Get the polynomial order of the elements.
   * @return ORDER
   */
  PROXY_HOST_DEVICE
  int getOrder() const final
  {
    return static_cast<int>(order_);  // Cast to int for consistency
  }

  /**
   * @brief Get the boundary type of a given node.
   * @param n Global node index
   * @return A combination of BoundaryFlag values
   */
  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(ScalarType n) const final
  {
    return static_cast<BoundaryFlag>(boundaries_t_[n]);
  }

  /**
   * @brief Compute the outward unit normal vector of an element face.
   * @param e Element index
   * @param dir Axis direction (0 = x, 1 = y, 2 = z)
   * @param face 1 = negative side, 2 = positive side in that direction
   * @param[out] v Output array (size 3) holding the normal vector
   */
  PROXY_HOST_DEVICE
  void faceNormal(ScalarType e, int dir, int face, FloatType v[3]) const final
  {
    // TODO: Implement actual face normal computation
    // throw std::runtime_error("FaceNormal not implemented");
    return;
  }

  /**
   * @brief Get the size of the domain in the specified dimension.
   * @param dim The dimension index (0 for X, 1 for Y, 2 for Z).
   * @return The size of the domain along the specified dimension.
   */
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

  /**
   * @brief Computes the minimum spacing between neighboring quadrature points.
   *
   * This function calculates the minimum Euclidean distance between adjacent
   * quadrature points in the spectral element mesh. Since all elements have
   * identical size and shape, the minimum spacing is computed by examining only
   * the first element (e=0), avoiding redundant calculations across all
   * elements.
   *
   * The algorithm checks spacing in three directions:
   * - i-direction: spacing between points (i, j, k) and (i+1, j, k)
   * - j-direction: spacing between points (i, j, k) and (i, j+1, k)
   * - k-direction: spacing between points (i, j, k) and (i, j, k+1)
   *
   * For each neighboring pair, the 3D Euclidean distance is computed:
   * distance = sqrt((x2-x1)² + (y2-y1)² + (z2-z1)²)
   *
   * @return The minimum spacing (in physical coordinates) between any two
   *         neighboring quadrature points in the mesh.
   *
   * @note This function assumes all elements are identical in size and shape.
   * @note Only checks direct neighbors along grid lines, not diagonal
   * neighbors.
   * @note Complexity: O(order³) instead of O(n_element × order³)
   */
  PROXY_HOST_DEVICE
  FloatType getMinSpacing() const final
  {
    FloatType minSpacing = std::numeric_limits<FloatType>::max();

    // Since all elements are the same size, only check the first element
    constexpr ScalarType e = 0;

    // Check i-direction spacing
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

    // Check j-direction spacing
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

    // Check k-direction spacing
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
    FloatType maxSpeedNode;
    FloatType maxSpeedElem;

    if (model_vp_node_.extent(0) > 0)
    {
      FIND_MAX(model_vp_node_, n_node_, maxSpeedNode);
    }
    else if (model_vp_element_.extent(0) > 0)
    {
      FIND_MAX(model_vp_element_, n_element_, maxSpeedElem);
    }
    else
    {
      throw std::runtime_error(
          "No model initialized (model unstruct getMaxSpeed).");
    }
    return max(maxSpeedElem, maxSpeedNode);
  }

 private:
  ScalarType order_;
  ScalarType n_element_;
  ScalarType n_node_;
  FloatType lx_, ly_, lz_;
  int n_points_per_element_;  // Added missing member
  bool isModelOnNodes_;

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

  VECTOR_REAL_VIEW boundaries_t_;
};
}  // namespace model

#endif  // SRC_MODEL_MODELAPI_INCLUDE_MODEL_UNSTRUCT_H_
