#ifndef BASE_MESH_
#define BASE_MESH_

#include <SEMmacros.hpp>
#include <dataType.hpp>

template <typename Coord, typename NodeIDX, typename ElementIDX, int ORDER>
class BaseMesh {

public:
  PROXY_HOST_DEVICE BaseMesh() {};

  PROXY_HOST_DEVICE
  BaseMesh(const ElementIDX &ex_in, const ElementIDX &ey_in,
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
  virtual ~BaseMesh() {}

  PROXY_HOST_DEVICE
  virtual Coord nodeCoordX(NodeIDX dofGlobal) const = 0;
  PROXY_HOST_DEVICE
  virtual Coord nodeCoordZ(NodeIDX dofGlobal) const = 0;
  PROXY_HOST_DEVICE
  virtual Coord nodeCoordY(NodeIDX dofGlobal) const = 0;

  PROXY_HOST_DEVICE
  virtual NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const = 0;

  PROXY_HOST_DEVICE
  ElementIDX getNumberOfElements() const { return ex * ey * ez; };

  PROXY_HOST_DEVICE
  NodeIDX getNumberOfNodes() const { return nx * ny * nz; };

  PROXY_HOST_DEVICE
  constexpr int getNumberOfPointsPerElement() const {
    return (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  }

  PROXY_HOST_DEVICE
  constexpr int getOrder() const { return ORDER; };

  int getDx() const { return hx; }
  int getDy() const { return hy; }
  int getDz() const { return hz; }

  ElementIDX getEx() const { return ex; }
  ElementIDX getEy() const { return ey; }
  ElementIDX getEz() const { return ez; }

  NodeIDX getNx() const { return nx; };
  NodeIDX getNy() const { return ny; };
  NodeIDX getNz() const { return nz; };

protected:
  ElementIDX ex, ey, ez; // Nb elements in each direction
  ElementIDX nx, ny, nz; // Nb nodes in each direction
  float lx, ly, lz;      // domain size
  float hx, hy, hz;      // element size
  int orderx, ordery, orderz, order;
};

#endif // BASE_MESH_
