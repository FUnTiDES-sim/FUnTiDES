#ifndef SEM_MESH_
#define SEM_MESH_

#include <SEMmacros.hpp>
#include <cmath>
#include <dataType.hpp>
#include <stdexcept>
using namespace std;

template <typename Coord, typename NodeIDX, typename ElementIDX> class SEMmesh {

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
    ny = ey * ordery + 1;
    nz = ez * orderz + 1;

    hx = lx / static_cast<float>(ex);
    hy = ly / static_cast<float>(ey);
    hz = lz / static_cast<float>(ez);
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordX(NodeIDX dofGlobal) const {
    int gx = dofGlobal % nx;
    return gx * hx / order;
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordZ(NodeIDX dofGlobal) const {
    int gz = (dofGlobal / nx) % nz;
    return gz * hz / order;
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordY(NodeIDX dofGlobal) const {
    int gy = dofGlobal / (nx * nz);
    return gy * hy / order;
  }

  PROXY_HOST_DEVICE
  NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const {
    ElementIDX elementI, elementJ, elementK;
    elementJ = e / (ex * ez);
    ElementIDX localIK = e - elementJ * ex * ez;
    elementI = localIK % ex;
    elementK = localIK / ex;

    ElementIDX elementOffset =
        elementI * order + elementK * order * nx + elementJ * order * nx * nz;

    NodeIDX dofGlobal = elementOffset + i + k * nx + j * nx * nz;
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

private:
  ElementIDX ex, ey, ez; // Nb elements in each direction
  ElementIDX nx, ny, nz; // Nb nodes in each direction
  float lx, ly, lz;      // domain size
  float hx, hy, hz;      // element size
  int orderx, ordery, orderz, order;
};
#endif // SEM_MESH_
