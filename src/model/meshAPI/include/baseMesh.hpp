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

  /**
  * @enum BoundaryFlag
  * @brief Flags representing the boundary condition type of a mesh node.
  *
  * This enumeration is used to mark nodes of the SEM mesh with specific
  * boundary properties. Multiple flags can be combined using bitwise OR.
  *
  * - `InteriorNode` : Node inside the computational domain.
  * - `Damping`      : Node in a damping boundary region.
  * - `Sponge`       : Node in a sponge layer region.
  * - `Surface`      : Node on a free surface boundary.
  * - `Ghost`        : Ghost node used for halo/exchange regions.
  *
  * @note Flags are implemented as bitfields, allowing combinations using `|`.
  *
  * Example:
  * @code
  * BoundaryFlag f = static_cast<BoundaryFlag>(Damping | Sponge);
  * @endcode
  */
  enum BoundaryFlag: uint8_t
  {
      InteriorNode = 0,
      Damping      = 1 << 0,
      Sponge       = 1 << 1,
      Surface      = 1 << 2,
      Ghost        = 1 << 3
  };

  /**
  * @brief Get the boundary type of a given node.
  *
  * This pure virtual function returns the boundary type of the specified node
  * in the SEM mesh.
  *
  * @param n Index of the node (`NodeIDX`) in the mesh.
  * @return A @ref BoundaryFlag indicating the boundary type (or `InteriorNode` if none).
  *
  * @note Must be implemented by derived classes for specific mesh/boundary setups.
  */
  PROXY_HOST_DEVICE
  virtual BoundaryFlag boundaryType(NodeIDX n) const = 0;

  /**
  * @brief Compute the outward normal vector of a given element face.
  *
  * For a given element `e`, direction `dir`, and face index `face`,
  * this function stores the **unit outward normal vector** of the face
  * into the array `v[3]`.
  *
  * - `dir` specifies the axis-aligned face group:
  *   - `0` → Front/Back faces (normal along ±x)
  *   - `1` → Left/Right faces (normal along ±y)
  *   - `2` → Bottom/Top faces (normal along ±z)
  *
  * - `face` specifies which side in the given `dir`:
  *   - `1` → The negative side (e.g., front, left, bottom)
  *   - `2` → The positive side (e.g., back, right, top)
  *
  * The computed vector `v` is of type `ModelType`, which is a floating-point
  * type (e.g., `float`, `double`).
  *
  * @param e Index of the element (`ElementIDX`) in the mesh.
  * @param dir Face direction axis (`0=x`, `1=y`, `2=z`).
  * @param face Which face in the given direction (`1=negative`, `2=positive`).
  * @param[out] v Output array of size 3 containing the face normal vector.
  *
  * @note Must be implemented by derived classes based on the mesh layout.
  *
  * @par Example:
  * @code
  * ModelType normal[3];
  * mesh->faceNormal(e, 0, 1, normal); // Get normal of the "front" face
  * @endcode
  */
  PROXY_HOST_DEVICE
  virtual void faceNormal(ElementIDX e, int dir, int face, ModelType v[3]) const = 0;
};

#endif // BASE_MESH_
