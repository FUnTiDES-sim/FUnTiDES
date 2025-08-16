#ifndef SEM_MESH_
#define SEM_MESH_

#include "baseMesh.hpp"
#include "CartesianParams.hpp"
#include "gllpoints.hpp"
#include <SEMmacros.hpp>
#include <cmath>
#include <dataType.hpp>
#include <stdexcept>

/**
 * @brief 3D Cartesian Spectral Element Method (SEM) mesh.
 *
 * This class implements a structured, Cartesian mesh using the Spectral Element
 * Method (SEM) in 3D space. It inherits from the templated `BaseMesh` class and
 * provides concrete implementations for coordinate mapping and global indexing.
 *
 * @tparam coord_t      Coordinate type (e.g., float or double)
 * @tparam index_t    Type used to index global nodes (e.g., int or std::size_t)
 * @tparam ORDER      Polynomial interpolation order per element
 *
 * @see BaseMesh
 */
template <typename coord_t, typename index_t, int ORDER>
class CartesianSEMmesh : public BaseMesh<coord_t, index_t> {
public:
  PROXY_HOST_DEVICE
  CartesianSEMmesh() {};

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
  CartesianSEMmesh (const typename BaseMesh<coord_t, index_t>::DataStruct & params)
  {
    auto const& p = dynamic_cast<const CartesianParams<coord_t, index_t>&>(params);

    order = p.order;

    ex = p.ex;
    ey = p.ey;
    ez = p.ez;

    lx = p.lx;
    ly = p.ly;
    lz = p.lz;


    nx = ex * ORDER + 1;
    ny = ey * ORDER + 1;
    nz = ez * ORDER + 1;

    hx = lx / static_cast<float>(ex);
    hy = ly / static_cast<float>(ey);
    hz = lz / static_cast<float>(ez);
  }

  PROXY_HOST_DEVICE ~CartesianSEMmesh(){};

  PROXY_HOST_DEVICE
  coord_t nodeCoord(index_t dofGlobal, int dim) const override final
  {
    // Calculate total number of nodes per dimension
    int nodesPerDim[3];
    nodesPerDim[0] = (ex * ORDER) + 1;  // nx total nodes in X
    nodesPerDim[1] = (ey * ORDER) + 1;  // ny total nodes in Y
    nodesPerDim[2] = (ez * ORDER) + 1;  // nz total nodes in Z

    // Convert global node index to 3D node indices (i, j, k)
    int k = dofGlobal / (nodesPerDim[0] * nodesPerDim[1]);
    int remainder = dofGlobal % (nodesPerDim[0] * nodesPerDim[1]);
    int j = remainder / nodesPerDim[0];
    int i = remainder % nodesPerDim[0];

    int nodeIdx[3] = {i, j, k};

    // Determine which element this node belongs to and local position within element
    int elemIdx = nodeIdx[dim] / ORDER;  // Element index in the requested dimension
    int localIdx = nodeIdx[dim] % ORDER; // Local node index within element (0 to ORDER)

    // Handle boundary case: if we're at the last node of an element (except the last element),
    // it's actually the first node of the next element
    if (localIdx == ORDER && elemIdx < (dim == 0 ? ex : (dim == 1 ? ey : ez)) - 1) {
        elemIdx++;
        localIdx = 0;
    }

    // Get the GLL point coordinate in reference element [-1, 1]
    coord_t gllPoint = GLLPoints<ORDER>::get(localIdx);

    // Map from reference element to physical element
    coord_t elementSize = (dim == 0) ? hx : ((dim == 1) ? hy : hz);
    coord_t elementStart = elemIdx * elementSize;

    // Transform from [-1, 1] to physical coordinates
    coord_t physicalCoord = elementStart + (gllPoint + 1.0) * elementSize * 0.5;

    return physicalCoord;
  }


  PROXY_HOST_DEVICE
  index_t globalNodeIndex(index_t e, int i, int j, int k) const override final
  {
    index_t elemZ = e / (ex * ey);
    index_t tmp   = e % (ex * ey);
    index_t elemY = tmp / ex;
    index_t elemX = tmp % ex;

    int ix = elemX * ORDER + i;
    int iy = elemY * ORDER + j;
    int iz = elemZ * ORDER + k;

    return ix + iy * nx + iz * nx * ny;
  }

  PROXY_HOST_DEVICE
  index_t getNumberOfElements() const override final
  {
    return this->ex * this->ey * this->ez;
  }

  PROXY_HOST_DEVICE
  coord_t getModelVpOnNodes(index_t n) const override final
  {
    return 1500;
  }

  PROXY_HOST_DEVICE
  coord_t getModelVpOnElement(index_t e) const override final
  {
    return 1500;
  }

  PROXY_HOST_DEVICE
  coord_t getModelRhoOnNodes(index_t n) const override final
  {
    return 1;
  }

  PROXY_HOST_DEVICE
  coord_t getModelRhoOnElement(index_t e) const override final
  {
    return 1;
  }


  PROXY_HOST_DEVICE
  index_t getNumberOfNodes() const override final
  {
    return this->nx * this->ny * this->nz;
  }

  PROXY_HOST_DEVICE
  int getNumberOfPointsPerElement() const override final
  {
    return (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  }

  PROXY_HOST_DEVICE
  int getOrder() const override final
  {
    return ORDER;
  }

  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(index_t n) const override final
  {
    return BoundaryFlag::InteriorNode;
  }

  PROXY_HOST_DEVICE
  void faceNormal(index_t e, int dir, int face, coord_t v[3]) const override final
  {
    return;
  }

  PROXY_HOST_DEVICE
  coord_t domainSize(int dim) const override final
  {
    if (dim == 0) return lx;
    if (dim == 1) return ly;
    return lz;
  }

  PROXY_HOST_DEVICE
  index_t elementFromCoordinate(coord_t x, coord_t y, coord_t z) const override final
  {
      // Calculate grid indices by scaling coordinates to grid space
      int i = static_cast<int>(x * ex / lx);
      int j = static_cast<int>(y * ey / ly);
      int k = static_cast<int>(z * ez / lz);

      // Calculate linear index using row-major ordering
      int index = i + ex * (j + ey * k);

      return index_t(index);
  }

#ifndef USE_KOKKOS
  VECTOR_REAL_VIEW
  extractXYSlice(const VECTOR_REAL_VIEW& array, index_t size, index_t z) const override final
  {
    int expected_size = size * size * size;

    // // Calculate slice parameters
    int slice_size = size * size;
    int start_index = z * slice_size;

    // // Create output view
    VECTOR_REAL_VIEW xy_slice = allocateVector<vectorReal>(slice_size);

    // Extract the slice using parallel_for
    // Kokkos::parallel_for("extract_xy_slice", slice_size, KOKKOS_LAMBDA(int i) {
    // LOOPHEAD(slice_size, i)
        // xy_slice(i) = array_1d(start_index + i);
    // LOOPEND
    // });

    return xy_slice;
  }

#else // USE_KOKKOS
  VECTOR_REAL_VIEW
  extractXYSlice(const VECTOR_REAL_VIEW& array, index_t size, index_t z)
  const override final
  {
      // Validate inputs
      if (z < 0 || z >= size) {
          Kokkos::abort("Z index out of bounds");
      }

      int slice_size = size * size;
      int start_index = z * slice_size;
      int end_index = start_index + slice_size;

      // Create subview (zero-copy operation)
      return Kokkos::subview(array, Kokkos::make_pair(start_index, end_index));
  }
#endif // USE_KOKKOS




private:
  index_t ex, ey, ez; // Nb elements in each direction
  index_t nx, ny, nz; // Nb nodes in each direction
  float lx, ly, lz;      // domain size
  float hx, hy, hz;      // element size
  int orderx, ordery, orderz, order;

};
#endif // SEM_MESH_
