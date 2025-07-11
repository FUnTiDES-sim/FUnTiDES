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
 * @tparam Coord      Type used to represent spatial coordinates ( real number)
 * @tparam NodeIDX    Type used for indexing global nodes (e.g., int or long)
 * @tparam ElementIDX Type used for indexing elements
 * @tparam ORDER      Polynomial order of the mesh elements
 */
template <typename Coord, typename NodeIDX, typename ElementIDX, int ORDER>
class BaseMesh {

public:
  /**
   * @brief Default constructor.
   *
   * Initializes an empty BaseMesh object.
   */
  PROXY_HOST_DEVICE BaseMesh() {};

  /**
   * @brief Constructs a structured mesh with element counts and physical sizes.
   *
   * @param ex_in   Number of elements in the X direction
   * @param ey_in   Number of elements in the Y direction
   * @param ez_in   Number of elements in the Z direction
   * @param lx_in   Physical length in the X direction
   * @param ly_in   Physical length in the Y direction
   * @param lz_in   Physical length in the Z direction
   * @param order   Polynomial interpolation order per element
   */
  PROXY_HOST_DEVICE
  BaseMesh(const ElementIDX &ex_in, const ElementIDX &ey_in,
           const ElementIDX &ez_in, const float &lx_in, const float &ly_in,
           const float &lz_in, const int order)
      : ex(ex_in), ey(ey_in), ez(ez_in), lx(lx_in), ly(ly_in), lz(lz_in),
        order(order), orderx(order), ordery(order), orderz(order) {
    nx = ex * orderx + 1;
    ny = ey * ordery + 1;
    nz = ez * orderz + 1;

    hx = lx / static_cast<float>(ex);
    hy = ly / static_cast<float>(ey);
    hz = lz / static_cast<float>(ez);
  }

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

  /**
   * @brief Get the total number of elements in the mesh.
   * @return ex * ey * ez
   */
  PROXY_HOST_DEVICE
  ElementIDX getNumberOfElements() const { return ex * ey * ez; };

  /**
   * @brief Get the total number of nodes in the mesh.
   * @return nx * ny * nz
   */
  PROXY_HOST_DEVICE
  NodeIDX getNumberOfNodes() const { return nx * ny * nz; };

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
   * @brief Get the physical spacing (hx) in the X direction.
   * @return Element width in X direction
   */
  int getDx() const { return hx; }

  /**
   * @brief Get the physical spacing (hy) in the Y direction.
   * @return Element width in Y direction
   */
  int getDy() const { return hy; }

  /**
   * @brief Get the physical spacing (hz) in the Z direction.
   * @return Element width in Z direction
   */
  int getDz() const { return hz; }

  /**
   * @brief Get the number of elements in the X direction.
   */
  ElementIDX getEx() const { return ex; }

  /**
   * @brief Get the number of elements in the Y direction.
   */
  ElementIDX getEy() const { return ey; }

  /**
   * @brief Get the number of elements in the Z direction.
   */
  ElementIDX getEz() const { return ez; }

  /**
   * @brief Get the number of nodes in the X direction.
   */
  NodeIDX getNx() const { return nx; };

  /**
   * @brief Get the number of nodes in the Y direction.
   */
  NodeIDX getNy() const { return ny; };

  /**
   * @brief Get the number of nodes in the Z direction.
   */
  NodeIDX getNz() const { return nz; };

protected:
  ElementIDX ex, ey, ez; // Nb elements in each direction
  ElementIDX nx, ny, nz; // Nb nodes in each direction
  float lx, ly, lz;      // domain size
  float hx, hy, hz;      // element size
  int orderx, ordery, orderz, order;
};

#endif // BASE_MESH_
