#ifndef SEM_MESH_
#define SEM_MESH_

#include "baseMesh.hpp"
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
 * @tparam ModelType  Type used for model parameters (e.g., float or double)
 * @tparam Coord      Coordinate type (e.g., float or double)
 * @tparam NodeIDX    Type used to index global nodes (e.g., int or std::size_t)
 * @tparam ElementIDX Type used to index elements
 * @tparam ORDER      Polynomial interpolation order per element
 *
 * @see BaseMesh
 */
template <typename ModelType, typename Coord, typename NodeIDX, typename ElementIDX, int ORDER>
class CartesianSEMmesh : public BaseMesh<ModelType, Coord, NodeIDX, ElementIDX, ORDER> {
public:
  PROXY_HOST_DEVICE
  CartesianSEMmesh() {};

  using Base = BaseMesh<ModelType, Coord, NodeIDX, ElementIDX, ORDER>;

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
  CartesianSEMmesh(const ElementIDX &ex_in, const ElementIDX &ey_in,
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

  PROXY_HOST_DEVICE ~CartesianSEMmesh(){};

  PROXY_HOST_DEVICE
  Coord nodeCoord(NodeIDX dofGlobal, int dim) const {
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
    Coord gllPoint = GLLPoints<ORDER>::get(localIdx);

    // Map from reference element to physical element
    Coord elementSize = (dim == 0) ? hx : ((dim == 1) ? hy : hz);
    Coord elementStart = elemIdx * elementSize;

    // Transform from [-1, 1] to physical coordinates
    Coord physicalCoord = elementStart + (gllPoint + 1.0) * elementSize * 0.5;

    return physicalCoord;
  }


  PROXY_HOST_DEVICE
  NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const {
    ElementIDX elemZ = e / (ex * ey);
    ElementIDX tmp   = e % (ex * ey);
    ElementIDX elemY = tmp / ex;
    ElementIDX elemX = tmp % ex;

    int ix = elemX * ORDER + i;
    int iy = elemY * ORDER + j;
    int iz = elemZ * ORDER + k;

    return ix + iy * nx + iz * nx * ny;
  }

  PROXY_HOST_DEVICE
  ElementIDX getNumberOfElements() const {
    return this->ex * this->ey * this->ez;
  }

  PROXY_HOST_DEVICE
  ModelType getModelVpOnNodes(NodeIDX n) const { return 1500; }

  PROXY_HOST_DEVICE
  ModelType getModelVpOnElement(ElementIDX e) const { return 1500; }

  PROXY_HOST_DEVICE
  ModelType getModelRhoOnNodes(NodeIDX n) const { return 1; }

  PROXY_HOST_DEVICE
  ModelType getModelRhoOnElement(ElementIDX e) const { return 1; }


  PROXY_HOST_DEVICE
  NodeIDX getNumberOfNodes() const { return this->nx * this->ny * this->nz; };

  PROXY_HOST_DEVICE
  constexpr int getNumberOfPointsPerElement() const {
    return (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  }

  PROXY_HOST_DEVICE
  constexpr int getOrder() const { return ORDER; }

  void extractXYslice(int k, vectorInt slice) const {
    int id = 0;
    for (int j = 0; j < this->ny; ++j) {
      for (int i = 0; i < this->nx; ++i) {
        NodeIDX nodeIdx = i + j * this->nx + k * this->nx * this->ny;
        slice(id) = nodeIdx;
      }
    }
  }

  PROXY_HOST_DEVICE
  typename Base::BoundaryFlag boundaryType(NodeIDX n) const
  {
    return Base::InteriorNode;
  }

  PROXY_HOST_DEVICE
  virtual void faceNormal(ElementIDX e, int dir, int face, ModelType v[3]) const {}

private:
  ElementIDX ex, ey, ez; // Nb elements in each direction
  ElementIDX nx, ny, nz; // Nb nodes in each direction
  float lx, ly, lz;      // domain size
  float hx, hy, hz;      // element size
  int orderx, ordery, orderz, order;

};
#endif // SEM_MESH_
