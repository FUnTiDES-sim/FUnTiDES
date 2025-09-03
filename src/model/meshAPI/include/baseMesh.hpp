#ifndef BASE_MESH_
#define BASE_MESH_

#include <SEMmacros.hpp>
#include <dataType.hpp>

/**
 * @enum BoundaryFlag
 * @brief Flags representing the boundary condition type of a mesh node.
 *
 * This enumeration is used to mark nodes of the SEM mesh with specific
 * boundary properties. Multiple flags can be combined using bitwise OR.
 */
enum BoundaryFlag: uint8_t
{
  InteriorNode = 0,          ///< Node inside the domain
  Damping      = 1 << 0,     ///< Node in damping boundary zone
  Sponge       = 1 << 1,     ///< Node in sponge layer
  Surface      = 1 << 2,     ///< Node on a free surface
  Ghost        = 1 << 3      ///< Ghost node for halo/exchange
};


/**
 * @brief Abstract base class representing a structured 3D mesh.
 *
 * This class defines the geometric and topological structure of a 3D mesh
 * using template parameters for coordinate types, node indices, element
 * indices, and polynomial order. It supports element-based and node-based
 * queries and defines an interface for accessing node coordinates and global
 * indices, as well as model parameters such as velocity and density.
 *
 * @tparam coord_t           Type used to represent spatial coordinates (real number)
 * @tparam index_t           Type used for indexing global nodes and elements
 */
template< typename coord_t, typename index_t >
class BaseMesh
{

public:
  /**
   * @brief Default constructor.
   *
   * Initializes an empty BaseMesh object.
   */
  PROXY_HOST_DEVICE BaseMesh() {}

  /**
   * @brief Virtual destructor.
   */
  PROXY_HOST_DEVICE ~BaseMesh() {}

  /**
   * @brief Get the coordinate of a global node in the given dimension.
   *
   * @param dofGlobal Global node index
   * @param dim Dimension index (0 = x, 1 = y, 2 = z)
   * @return Coordinate value in the specified dimension
   */
  PROXY_HOST_DEVICE
  virtual coord_t nodeCoord( index_t dofGlobal, int dim ) const = 0;

  /**
   * @brief Get the global node index for a local element-node triplet.
   *
   * @param e Element index
   * @param i Local i-index in the element
   * @param j Local j-index in the element
   * @param k Local k-index in the element
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  virtual index_t globalNodeIndex( index_t e, int i, int j, int k ) const = 0;

  /**
   * @brief Get the P-wave velocity value at a global node.
   *
   * @param n Global node index
   * @return Model P-wave velocity value at the node
   */
  PROXY_HOST_DEVICE
  virtual coord_t getModelVpOnNodes( index_t n ) const = 0;

  /**
   * @brief Get the average P-wave velocity value on a given element.
   *
   * @param e Element index
   * @return Model P-wave velocity value for the element
   */
  PROXY_HOST_DEVICE
  virtual coord_t getModelVpOnElement( index_t e ) const = 0;

  /**
   * @brief Get the density value at a global node.
   *
   * @param n Global node index
   * @return Model density value at the node
   */
  PROXY_HOST_DEVICE
  virtual coord_t getModelRhoOnNodes( index_t n ) const = 0;

  /**
   * @brief Get the average density value on a given element.
   *
   * @param e Element index
   * @return Model density value for the element
   */
  PROXY_HOST_DEVICE
  virtual coord_t getModelRhoOnElement( index_t e ) const = 0;

  /**
   * @brief Get the total number of elements in the mesh.
   * @return Total element count (typically ex * ey * ez)
   */
  PROXY_HOST_DEVICE
  virtual index_t getNumberOfElements() const = 0;

  /**
   * @brief Get the total number of global nodes in the mesh.
   * @return Total node count (typically nx * ny * nz)
   */
  PROXY_HOST_DEVICE
  virtual index_t getNumberOfNodes() const = 0;

  /**
   * @brief Get the number of interpolation points per element.
   *
   * This is equal to (ORDER + 1)^3 due to tensor-product basis in SEM.
   *
   * @return Number of interpolation points in one element
   */
  PROXY_HOST_DEVICE
  virtual int getNumberOfPointsPerElement() const = 0;

  /**
   * @brief Get the polynomial order of the elements.
   *
   * @return ORDER (template parameter)
   */
  PROXY_HOST_DEVICE
  virtual int getOrder() const = 0;

  /**
   * @brief Get the boundary type of a given node.
   *
   * @param n Global node index
   * @return A combination of BoundaryFlag values
   */
  PROXY_HOST_DEVICE
  virtual BoundaryFlag boundaryType( index_t n ) const = 0;

  /**
   * @brief Compute the outward unit normal vector of an element face.
   *
   * @param e Element index
   * @param dir Axis direction (0 = x, 1 = y, 2 = z)
   * @param face 1 = negative side, 2 = positive side in that direction
   * @param[out] v Output array (size 3) holding the normal vector
   */
  PROXY_HOST_DEVICE
  virtual void faceNormal( index_t e, int dir, int face, coord_t v[3] ) const = 0;

  /**
   * @brief Get the size of the domain in the specified dimension.
   *
   * This function returns the size of the domain along the given dimension.
   * Dimensions are typically 0 (X), 1 (Y), or 2 (Z).
   *
   * @param dim The dimension index (0 for X, 1 for Y, 2 for Z).
   * @return The size of the domain along the specified dimension.
   */
  PROXY_HOST_DEVICE
  virtual coord_t domainSize( int dim ) const = 0;

  /**
   * @brief Get the element index at the given coordinate.
   *
   * This function returns the index of the element located at the
   * specified (x, y, z) coordinates within the domain.
   *
   * @param x The X coordinate.
   * @param y The Y coordinate.
   * @param z The Z coordinate.
   * @return The element index corresponding to the given coordinates.
   */
  PROXY_HOST_DEVICE
  virtual index_t elementFromCoordinate( coord_t x, coord_t y, coord_t z ) const = 0;

  /**
   * Extract an XY slice from a Kokkos 1D View representing a 3D cubic array
   *
   * For Kokkos:
   * Extract XY slice using subview (more efficient, zero-copy)
   * Returns a subview that shares memory with the original array
   *
   * @param array_1d: Input Kokkos 1D View containing 3D data (stored in X-Y-Z order)
   * @param size: Size of each dimension (assuming cubic array: size x size x size)
   * @param z: Z-level to extract (0 to size-1)
   * @return: Kokkos 1D View containing the XY slice
   */
  virtual
  VECTOR_REAL_VIEW
  extractXYSlice( const VECTOR_REAL_VIEW & array, index_t size, index_t z ) const
  {
    std::cout << "ExctractXYSlice is not implemented in MeshBase" << std::endl;
    throw std::runtime_error( "ExtractXYSlice not implemented" );  // Proper exception
  };

};

#endif // BASE_MESH_
