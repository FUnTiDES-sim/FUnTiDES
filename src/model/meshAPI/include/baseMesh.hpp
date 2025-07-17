#ifndef BASE_MESH_
#define BASE_MESH_

#include <SEMmacros.hpp>
#include <dataType.hpp>

/**
 * @brief Abstract base class representing a structured 3D mesh.
 *
 * This class defines the geometric and topological structure of a 3D mesh
 * using template parameters for coordinate types, node indices, element
 * indices, and polynomial order. It supports element-based and node-based
 * queries and defines an interface for accessing node coordinates and global
 * indices.
 *
 * @tparam ModelType  Type used for model parameters (e.g., float or double)
 * @tparam Coord      Type used to represent spatial coordinates ( real number)
 * @tparam NodeIDX    Type used for indexing global nodes (e.g., int or long)
 * @tparam ElementIDX Type used for indexing elements
 * @tparam ORDER      Polynomial order of the mesh elements
 */
template <typename ModelType, typename Coord, typename NodeIDX, typename ElementIDX, int ORDER>
class BaseMesh {

public:
  /**
   * @brief Default constructor.
   *
   * Initializes an empty BaseMesh object.
   */
  PROXY_HOST_DEVICE BaseMesh() {};

  /**
   * @brief Virtual destructor.
   */
  PROXY_HOST_DEVICE
  virtual ~BaseMesh() {}
  /**
   * @brief Get the X coordinate of a global node.
   * @param dofGlobal Global node index
   * @return X coordinate of the node
   */
  PROXY_HOST_DEVICE
  virtual Coord nodeCoordX(NodeIDX dofGlobal) const = 0;

  /**
   * @brief Get the Y coordinate of a global node.
   * @param dofGlobal Global node index
   * @return Y coordinate of the node
   */
  PROXY_HOST_DEVICE
  virtual Coord nodeCoordY(NodeIDX dofGlobal) const = 0;

  /**
   * @brief Get the Z coordinate of a global node.
   * @param dofGlobal Global node index
   * @return Z coordinate of the node
   */
  PROXY_HOST_DEVICE
  virtual Coord nodeCoordZ(NodeIDX dofGlobal) const = 0;

  /**
   * @brief Get the global node index for a local element-node triplet.
   *
   * @param e Element index
   * @param i Local i-index in the element (0 ≤ i ≤ ORDER)
   * @param j Local j-index in the element (0 ≤ j ≤ ORDER)
   * @param k Local k-index in the element (0 ≤ k ≤ ORDER)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  virtual NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const = 0;

  PROXY_HOST_DEVICE
  virtual ModelType getModelVpOnNodes(NodeIDX n) const = 0;

  PROXY_HOST_DEVICE
  virtual ModelType getModelVpOnElement(ElementIDX e) const = 0;

  PROXY_HOST_DEVICE
  virtual ModelType getModelRhoOnNodes(NodeIDX n) const = 0;

  PROXY_HOST_DEVICE
  virtual ModelType getModelRhoOnElement(ElementIDX e) const = 0;

  /**
   * @brief Get the total number of elements in the mesh.
   * @return ex * ey * ez
   */
  PROXY_HOST_DEVICE
  virtual ElementIDX getNumberOfElements() const = 0;

  /**
   * @brief Get the total number of nodes in the mesh.
   * @return nx * ny * nz
   */
  PROXY_HOST_DEVICE
  virtual NodeIDX getNumberOfNodes() const = 0;

  /**
   * @brief Get the number of interpolation points per element.
   * @return (ORDER + 1)^3
   */
  PROXY_HOST_DEVICE
  constexpr int getNumberOfPointsPerElement() const {
    return (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  }

  /**
   * @brief Get the polynomial order of the elements.
   * @return ORDER (template constant)
   */
  PROXY_HOST_DEVICE
  constexpr int getOrder() const { return ORDER; };

};

#endif // BASE_MESH_
