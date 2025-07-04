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
          const float &lz_in, const int &order_in)
      : ex(ex_in), ey(ey_in), ez(ez_in), lx(lx_in), ly(ly_in), lz(lz_in),
        order(order_in), orderx(order_in), ordery(order_in), orderz(order_in) {
    nx = ex * orderx + 1;
    ny = (ey == 0) ? 1 : ey * ordery + 1;
    nz = ez * orderz + 1;

    hx = lx / static_cast<float>(ex);
    hy = (ey == 0) ? 1.0f : ly / static_cast<float>(ey);
    hz = lz / static_cast<float>(ez);
  }

  PROXY_HOST_DEVICE
  Coord nodeCoordX(NodeIDX n);
  PROXY_HOST_DEVICE
  Coord nodeCoordY(NodeIDX n);
  PROXY_HOST_DEVICE
  Coord nodeCoordZ(NodeIDX n);

  PROXY_HOST_DEVICE
  NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const;

  PROXY_HOST_DEVICE
  ElementIDX getNumberOfElements() const { return ex * ey * ez; };

  PROXY_HOST_DEVICE
  NodeIDX getNumberOfNodes() const { return nx * ny * nz; };

  PROXY_HOST_DEVICE
  int getNumberOfPointsPerElement() const {
    return (order + 1) * (order + 1) * (order + 1);
  }

  PROXY_HOST_DEVICE
  int getOrder() const { return order; };

  PROXY_HOST_DEVICE
  int getModel(ElementIDX e) const { return 1500; };

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
