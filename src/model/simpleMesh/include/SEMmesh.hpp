#ifndef SEM_MESH_
#define SEM_MESH_

#include <SEMmacros.hpp>
#include <cmath>
#include <dataType.hpp>
#include <stdexcept>
using namespace std;

template <typename Coord, typename NodeIDX, typename ElementIDX> class SEMmesh {
private:
  ElementIDX ex, ey, ez;
  ElementIDX nx, ny, nz;
  float lx, ly, lz;
  float hx, hy, hz;
  int orderx, ordery, orderz, order;

public:
  PROXY_HOST_DEVICE SEMmesh() {};
  PROXY_HOST_DEVICE ~SEMmesh(){};

  PROXY_HOST_DEVICE
  SEMmesh(const ElementIDX &ex_in, const ElementIDX &ey_in,
          const ElementIDX &ez_in, const float &lx_in, const float &ly_in,
          const float &lz_in, const int order)
      : ex(ex_in), ey(ey_in), ez(ez_in), lx(lx_in), ly(ly_in), lz(lz_in),
        order(order), orderx(order), ordery(order), orderz(order) {
    nx = ex * orderx + 1;
    ny = (ey == 0) ? 1 : ey * ordery + 1;
    nz = ez * orderz + 1;

    hx = lx / static_cast<float>(ex);
    hy = (ey == 0) ? 1.0f : ly / static_cast<float>(ey);
    hz = lz / static_cast<float>(ez);
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordX(NodeIDX n) const { // get element n0 coordinate
    ElementIDX e = n / ((order + 1) * (order + 1) * (order + 1));
    int elementX = e % ex;
    Coord n0 = elementX * hx;
    // get node coordinate offset
    int dofIDX = n % ((order + 1) * (order + 1) * (order + 1));
    int dofx = dofIDX % (order + 1);
    int nodeXOffset =
        n0 + (hx / (order + 1) * dofx); // TODO: Takes GLL coordinate instead

    return nodeXOffset + n0;
  };
  PROXY_HOST_DEVICE
  Coord nodeCoordY(NodeIDX n) const { // Total number of DoFs per element
    int dofsPerElem = (order + 1) * (order + 1) * (order + 1);

    // Identify the element index this node belongs to
    ElementIDX e = n / dofsPerElem;

    // Compute the Y index of the element in the global mesh
    int elementY =
        (e / ex) % ey; // Assumes row-major ordering (X fastest, then Y, then Z)

    // Compute the origin coordinate of the element in Y direction
    Coord n0 = elementY * hy;

    // Determine the local node index within the element
    int dofIDX = n % dofsPerElem;

    // Extract the Y DoF index from local index
    int dofy = (dofIDX / (order + 1)) % (order + 1);

    // Compute node offset inside the element in Y direction
    int nodeYOffset =
        hy / (order + 1) * dofy; // TODO: Use GLL points if required

    return nodeYOffset + n0;
  };
  PROXY_HOST_DEVICE
  Coord nodeCoordZ(NodeIDX n) const { // get element n0 coordinate
    ElementIDX e = n / ((order + 1) * (order + 1) * (order + 1));
    int elementZ = e / ez;
    Coord n0 = elementZ * hz;
    // get node coordinate offset
    int dofIDX = n % ((order + 1) * (order + 1) * (order + 1));
    int dofz = dofIDX / (order + 1);
    int nodeZOffset =
        n0 + (hz / (order + 1) * dofz); // TODO: Takes GLL coordinate instead

    return nodeZOffset + n0;
  };

  PROXY_HOST_DEVICE
  NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const {
    ElementIDX elementI, elementJ, elementK;
    elementJ = e / (ex * ez);
    ElementIDX localIK = e - elementJ * ex * ez;
    elementI = localIK % ex;
    elementK = localIK / ex;

    ElementIDX elementOffset =
        elementI * order + elementK * order * nx + elementJ * order * nx * nz;

    NodeIDX dofGlobal = elementOffset + i + j * nx + k * nx * nz;
    return dofGlobal;
  };

  PROXY_HOST_DEVICE
  ElementIDX getNumberOfElements() const { return ex * ey * ez; };

  PROXY_HOST_DEVICE
  NodeIDX getNumberOfNodes() const { return nx * ny * nz; };

  PROXY_HOST_DEVICE
  constexpr int getNumberOfPointsPerElement() const {
    // return (order + 1) * (order + 1) * (order + 1);
    return 27; // TODO: This is temporary hack
  }

  PROXY_HOST_DEVICE
  constexpr int getOrder() const { return 2; };
  // TODO X Y Z
  // TODO Change constexpr

  PROXY_HOST_DEVICE
  int getModel(ElementIDX e) const { return 1500; };

  void extractXYslice(int k, vectorInt slice) const {
    int id = 0;
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        NodeIDX nodeIdx = i + j * nx + k * nx * ny;
        slice(id) = nodeIdx;
      }
    }
  }

  int getDx() const { return hx; }
  int getDy() const { return hy; }
  int getDz() const { return hz; }

  ElementIDX getEx() const { return ex; }
  ElementIDX getEy() const { return ey; }
  ElementIDX getEz() const { return ez; }

  NodeIDX getNx() const { return nx; };
  NodeIDX getNy() const { return ny; };
  NodeIDX getNz() const { return nz; };
};
#endif // SEM_MESH_
