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
  Coord nodeCoordX(NodeIDX dofGlobal) const {
    // int gx = dofGlobal % this->nx;
    // return gx * this->hx / this->order;
    constexpr int nodesPerElem = ORDER + 1;
    int elemId = dofGlobal / nodesPerElem;
    int localId = dofGlobal % nodesPerElem;

    double elemStart = elemId * hx;
    double elemEnd   = (elemId + 1) * hx;

    // Map reference [-1,1] -> [elemStart, elemEnd]
    double xi = GLLPoints<ORDER>::get(localId);

    // Affine map from [-1,1] to [elemStart, elemEnd]
    auto physicalX = 0.5 * ( (1 - xi) * elemStart + (1 + xi) * elemEnd );
    return physicalX;
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordZ(NodeIDX dofGlobal) const {
    int gz = (dofGlobal / this->nx) % this->nz;
    return gz * this->hz / this->order;
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordY(NodeIDX dofGlobal) const {
    int gy = dofGlobal / (this->nx * this->nz);
    return gy * this->hy / this->order;
  }

  PROXY_HOST_DEVICE
  NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const {
    ElementIDX elementI, elementJ, elementK;
    elementJ = e / (this->ex * this->ez);
    ElementIDX localIK = e - elementJ * this->ex * this->ez;
    elementI = localIK % this->ex;
    elementK = localIK / this->ex;

    ElementIDX elementOffset = elementI * this->order +
                               elementK * this->order * this->nx +
                               elementJ * this->order * this->nx * this->nz;

    NodeIDX dofGlobal =
        elementOffset + i + k * this->nx + j * this->nx * this->nz;
    return dofGlobal;
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
