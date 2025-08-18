#ifndef SRC_MODEL_CARTESIANSEMMESH_INCLUDE_CARTESIANSEMMESH_HPP_
#define SRC_MODEL_CARTESIANSEMMESH_INCLUDE_CARTESIANSEMMESH_HPP_

#include <cmath>
#include <stdexcept>
#include <dataType.hpp>
#include <SEMmacros.hpp>
#include "../../mesh_api/include/baseMesh.hpp"
#include "./cartesian_params.h"
#include "./gllpoints.hpp"

/**
 * @brief 3D Cartesian Spectral Element Method (SEM) mesh.
 *
 * This class implements a structured, Cartesian mesh using the Spectral Element
 * Method (SEM) in 3D space. It inherits from the templated `BaseMesh` class and
 * provides concrete implementations for coordinate mapping and global indexing.
 *
 * @tparam Coord      Coordinate type (e.g., float or double)
 * @tparam Index    Type used to index global nodes (e.g., int or std::size_t)
 * @tparam ORDER      Polynomial interpolation order per element
 *
 * @see BaseMesh
 */
template <typename Coord, typename Index, int ORDER>
class CartesianSEMmesh : public mesh_base::BaseMesh<Coord, Index> {
 public:
  PROXY_HOST_DEVICE
  CartesianSEMmesh() {}

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
  explicit CartesianSEMmesh(
      const typename mesh_base::BaseMesh<Coord, Index>::DataStruct& params) {
    auto const& p =
        dynamic_cast<const CartesianParams<Coord, Index>&>(params);

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

  PROXY_HOST_DEVICE ~CartesianSEMmesh() {}

  PROXY_HOST_DEVICE
  Coord nodeCoord(Index dofGlobal, int dim) const final {
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

    int elemIdx = nodeIdx[dim] / ORDER;
    int localIdx = nodeIdx[dim] % ORDER;

    if (localIdx == ORDER
        && elemIdx < (dim == 0 ? ex : (dim == 1 ? ey : ez)) - 1) {
      elemIdx++;
      localIdx = 0;
    }

    // Get the GLL point coordinate in reference element [-1, 1]
    Coord gllPoint = GLLPoints<ORDER>::get(localIdx);

    // Map from reference element to physical element
    Coord elementSize = (dim == 0) ? hx : ((dim == 1) ? hy : hz);
    Coord elementStart = elemIdx * elementSize;

    // Transform from [-1, 1] to physical coordinates
    Coord physicalCoord = elementStart + (gllPoint + 1.0) * elementSize * 0.5;

    return physicalCoord;
  }


  PROXY_HOST_DEVICE
  Index globalNodeIndex(Index e, int i, int j, int k) const final {
    Index elemZ = e / (ex * ey);
    Index tmp   = e % (ex * ey);
    Index elemY = tmp / ex;
    Index elemX = tmp % ex;

    int ix = elemX * ORDER + i;
    int iy = elemY * ORDER + j;
    int iz = elemZ * ORDER + k;

    return ix + iy * nx + iz * nx * ny;
  }

  PROXY_HOST_DEVICE
  Index getNumberOfElements() const final {
    return this->ex * this->ey * this->ez;
  }

  PROXY_HOST_DEVICE
  Coord getModelVpOnNodes(Index n) const final {
    return 1500;
  }

  PROXY_HOST_DEVICE
  Coord getModelVpOnElement(Index e) const final {
    return 1500;
  }

  PROXY_HOST_DEVICE
  Coord getModelRhoOnNodes(Index n) const final {
    return 1;
  }

  PROXY_HOST_DEVICE
  Coord getModelRhoOnElement(Index e) const final {
    return 1;
  }


  PROXY_HOST_DEVICE
  Index getNumberOfNodes() const final {
    return this->nx * this->ny * this->nz;
  }

  PROXY_HOST_DEVICE
  int getNumberOfPointsPerElement() const final {
    return (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  }

  PROXY_HOST_DEVICE
  int getOrder() const final {
    return ORDER;
  }

  PROXY_HOST_DEVICE
  mesh_base::BoundaryFlag boundaryType(Index n) const final {
    return mesh_base::BoundaryFlag::InteriorNode;
  }

  PROXY_HOST_DEVICE
  void faceNormal(Index e, int dir,
                  int face, Coord v[3]) const final {
    return;
  }

  PROXY_HOST_DEVICE
  Coord domainSize(int dim) const final {
    if (dim == 0) return lx;
    if (dim == 1) return ly;
    return lz;
  }

  PROXY_HOST_DEVICE
  Index elementFromCoordinate(Coord x, Coord y, Coord z) const final {
      // Calculate grid indices by scaling coordinates to grid space
      int i = static_cast<int>(x * ex / lx);
      int j = static_cast<int>(y * ey / ly);
      int k = static_cast<int>(z * ez / lz);

      // Calculate linear index using row-major ordering
      int index = i + ex * (j + ey * k);

      return Index(index);
  }

#ifndef USE_KOKKOS
  VECTOR_REAL_VIEW
  extractXYSlice(const VECTOR_REAL_VIEW& array,
                 Index size, Index z) const final {
    int expected_size = size * size * size;

    // Calculate slice parameters
    int slice_size = size * size;
    int start_index = z * slice_size;

    // Create output view
    VECTOR_REAL_VIEW xy_slice = allocateVector<vectorReal>(slice_size);

    return xy_slice;
  }

#else  // USE_KOKKOS
  VECTOR_REAL_VIEW
  extractXYSlice(const VECTOR_REAL_VIEW& array, Index size, Index z)
  const final {
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
#endif  // USE_KOKKOS

 private:
  Index ex, ey, ez;  // Nb elements in each direction
  Index nx, ny, nz;  // Nb nodes in each direction
  float lx, ly, lz;      // domain size
  float hx, hy, hz;      // element size
  int orderx, ordery, orderz, order;
};
#endif  // SRC_MODEL_CARTESIANSEMMESH_INCLUDE_CARTESIANSEMMESH_HPP_
