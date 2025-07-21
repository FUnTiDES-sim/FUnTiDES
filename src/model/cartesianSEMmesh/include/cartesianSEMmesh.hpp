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
       // 1) nb noeuds par élément
    int nodesPerElem = (order + 1) * (order + 1) * (order + 1);

    // 2) décomposition en élément + noeud local
    int elemIndex  = dofGlobal / nodesPerElem;
    int localIndex = dofGlobal % nodesPerElem;

    // --- élément l,m,n (même ordre que l’original) ---
    int l = elemIndex % ex;             // élément en X
    int m = (elemIndex / ex) % ey;      // élément en Y
    int n = elemIndex / (ex * ey);      // élément en Z

    // --- noeud local i,j,k ---
    int i = localIndex % (order + 1);
    int j = (localIndex / (order + 1)) % (order + 1);
    int k = localIndex / ((order + 1) * (order + 1));

    // choisir l’élément + noeud local + pas sur l’axe voulu
    int elemAxis  = 0;
    int localAxis = 0;
    float hAxis   = 0.0f;

    if (dim == 0) { // X
        elemAxis  = l;
        localAxis = i;
        hAxis     = hx;
    } else if (dim == 1) { // Y
        elemAxis  = m;
        localAxis = j;
        hAxis     = hy;
    } else { // Z
        elemAxis  = n;
        localAxis = k;
        hAxis     = hz;
    }
    float xi = GLLPoints<ORDER>::get(localAxis);
    float x0 = elemAxis * hAxis;
    float x1 = (elemAxis + 1) * hAxis;

    float b = 0.5f * (x0 + x1);  // centre
    float a = 0.5f * (x1 - x0);  // demi-longueur = h/2

    return a * xi + b;

  }


  PROXY_HOST_DEVICE
  NodeIDX globalNodeIndex(ElementIDX e, int i, int j, int k) const {
    ElementIDX elementI, elementJ, elementK;
    elementK = e / (this->ex * this->ez);
    ElementIDX localIJ = e - elementK * this->ex * this->ez;
    elementI = localIJ % this->ex;
    elementJ = localIJ / this->ex;

    ElementIDX elementOffset = elementI * this->order +
                               elementJ * this->order * this->nx +
                               elementK * this->order * this->nx * this->nz;

    NodeIDX dofGlobal =
        elementOffset + i + j * this->nx + k * this->nx * this->nz;
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
