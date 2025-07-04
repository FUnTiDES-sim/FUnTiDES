#include "SEMmesh.hpp"

template <typename Coord, typename NodeIDX, typename ElementIDX>
PROXY_HOST_DEVICE Coord
SEMmesh<Coord, NodeIDX, ElementIDX>::nodeCoordX(NodeIDX n) {
  // get element n0 coordinate
  ElementIDX e = n / ((order + 1) * (order + 1) * (order + 1));
  int elementX = e % ex;
  Coord n0 = elementX * hx;
  // get node coordinate offset
  int dofIDX = n % ((order + 1) * (order + 1) * (order + 1));
  int dofx = dofIDX % (order + 1);
  int nodeXOffset =
      n0 + (hx / (order + 1) * dofx); // TODO: Takes GLL coordinate instead

  return nodeXOffset + n0;
}

template <typename Coord, typename NodeIDX, typename ElementIDX>
PROXY_HOST_DEVICE Coord
SEMmesh<Coord, NodeIDX, ElementIDX>::nodeCoordZ(NodeIDX n) {
  // get element n0 coordinate
  ElementIDX e = n / ((order + 1) * (order + 1) * (order + 1));
  int elementZ = e / ez;
  Coord n0 = elementZ * hz;
  // get node coordinate offset
  int dofIDX = n % ((order + 1) * (order + 1) * (order + 1));
  int dofz = dofIDX / (order + 1);
  int nodeZOffset =
      n0 + (hz / (order + 1) * dofz); // TODO: Takes GLL coordinate instead

  return nodeZOffset + n0;
}

template <typename Coord, typename NodeIDX, typename ElementIDX>
PROXY_HOST_DEVICE Coord
SEMmesh<Coord, NodeIDX, ElementIDX>::nodeCoordY(NodeIDX n) {
  return 0;
}

template <typename Coord, typename NodeIDX, typename ElementIDX>
PROXY_HOST_DEVICE NodeIDX SEMmesh<Coord, NodeIDX, ElementIDX>::globalNodeIndex(
    ElementIDX e, int i, int j, int k) const {
  ElementIDX elementI, elementJ, elementK;
  elementJ = e / (ex * ez);
  ElementIDX localIK = e - elementJ * ex * ez;
  elementI = localIK % ex;
  elementK = localIK / ex;

  ElementIDX elementOffset =
      elementI * order + elementK * order * nx + elementJ * order * nx * nz;

  NodeIDX dofGlobal = elementOffset + i + j * nx + k * nx * nz;
  return dofGlobal;
}

template class SEMmesh<float, int, int>;
