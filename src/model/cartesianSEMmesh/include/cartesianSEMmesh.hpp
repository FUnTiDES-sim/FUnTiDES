#ifndef SEM_MESH_
#define SEM_MESH_

#include "baseMesh.hpp"
#include <SEMmacros.hpp>
#include <cmath>
#include <dataType.hpp>

template <typename Coord, typename NodeIDX, typename ElementIDX, int ORDER>
class CartesianSEMmesh : public BaseMesh<Coord, NodeIDX, ElementIDX, ORDER> {
public:
  PROXY_HOST_DEVICE CartesianSEMmesh() {};
  PROXY_HOST_DEVICE ~CartesianSEMmesh(){};

  PROXY_HOST_DEVICE
  CartesianSEMmesh(const ElementIDX &ex_in, const ElementIDX &ey_in,
                   const ElementIDX &ez_in, const float &lx_in,
                   const float &ly_in, const float &lz_in, const int order)
      : BaseMesh<Coord, NodeIDX, ElementIDX, ORDER>(ex_in, ey_in, ez_in, lx_in,
                                                    ly_in, lz_in, order) {}

  PROXY_HOST_DEVICE
  Coord nodeCoordX(NodeIDX dofGlobal) const {
    int gx = dofGlobal % this->nx;
    return gx * this->hx / this->order;
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
  NodeIDX getNumberOfNodes() const { return this->nx * this->ny * this->nz; };

  PROXY_HOST_DEVICE
  constexpr int getNumberOfPointsPerElement() const {
    return (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  }

  PROXY_HOST_DEVICE
  constexpr int getOrder() const { return ORDER; }
  // TODO X Y Z

  PROXY_HOST_DEVICE
  int getModel(ElementIDX e) const { return 1500; };

  void extractXYslice(int k, vectorInt slice) const {
    int id = 0;
    for (int j = 0; j < this->ny; ++j) {
      for (int i = 0; i < this->nx; ++i) {
        NodeIDX nodeIdx = i + j * this->nx + k * this->nx * this->ny;
        slice(id) = nodeIdx;
      }
    }
  }
};
#endif // SEM_MESH_
