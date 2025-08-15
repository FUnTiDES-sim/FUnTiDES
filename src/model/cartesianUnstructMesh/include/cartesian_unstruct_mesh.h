#ifndef SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
#define SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_

#include "baseMesh.hpp"
#include <dataType.hpp>
#include "commonMacros.hpp"

template <typename Coord, typename Index, int order>
class CartesianUnstructMesh : public BaseMesh<Coord, Index> {
 public:
  PROXY_HOST_DEVICE
  CartesianUnstructMesh() { }

  PROXY_HOST_DEVICE
    ~CartesianUnstructMesh() { }

  PROXY_HOST_DEVICE
  Coord nodeCoord(Index dofGlobalIndex, int dim) const final {
    switch (dim) {
      case 0: {
        return node_coord_x[dofGlobalIndex];
      }
      case 1: {
        return node_coord_y[dofGlobalIndex];
      }
      case 2: {
        return node_coord_z[dofGlobalIndex];
      }
      default:
        return -1;
    }
  }

  PROXY_HOST_DEVICE
  Index globalNodeIndex(Index elementIndex,
                        int i,
                        int j,
                        int k) const final {
    // TODO: Implem
    return -1;
  }

  PROXY_HOST_DEVICE
  Index getNumberOfElements() const final {
    return this->ex_ * this->ey_ * this->ez_;
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
    Index nx = this->ex_ * order + 1;
    Index ny = this->ey_ * order + 1;
    Index nz = this->ez_ * order + 1;
    return nx * ny * nz;
  }

  PROXY_HOST_DEVICE
  int getNumberOfPointsPerElement() const final {
    return (order + 1) * (order + 1) * (order + 1);
  }

  PROXY_HOST_DEVICE
  int getOrder() const final {
    return order;
  }

  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(Index n) const final {
    return BoundaryFlag::InteriorNode;
  }

  PROXY_HOST_DEVICE
  void faceNormal(Index e,
                  int dir,
                  int face,
                  Coord v[3]) const final {
    return;
  }

  PROXY_HOST_DEVICE
  Coord domainSize(int dim) const final {
    switch (dim) {
      case 0 :
        return lx_;
      case 1:
        return ly_;
      case 2:
        return lz_;
      default:
        return -1;
    }
  }

  PROXY_HOST_DEVICE
  Index elementFromCoordinate(Coord x, Coord y, Coord z) const final {
      int i = static_cast<int>(x * ex_ / lx_);
      int j = static_cast<int>(y * ey_ / ly_);
      int k = static_cast<int>(z * ez_ / lz_);

      // Calculate linear index using row-major ordering
      auto index = i + ex_ * (j + ey_ * k);

      return index;
  }

  VECTOR_REAL_VIEW
  extractXYSlice(const VECTOR_REAL_VIEW& array,
                 Index size,
                 Index z) const final {
#ifndef USE_KOKKOS
      throw "ExtractXYSlice is not working "
            "without kokkos (Cartesian Unstruct Mesh)"
#endif  // USE_KOKKOS
      if (z < 0 || z >= size) {
          Kokkos::abort("Z index out of bounds");
      }

      int slice_size = size * size;
      int start_index = z * slice_size;
      int end_index = start_index + slice_size;

      // Create subview (zero-copy operation)
      return Kokkos::subview(array, Kokkos::make_pair(start_index, end_index));
  }

 private:
  Index ex_, ey_, ez_;
  Coord lx_, ly_, lz_;

  ARRAY_REAL_VIEW global_node_index_;
  VECTOR_REAL_VIEW node_coord_x;
  VECTOR_REAL_VIEW node_coord_y;
  VECTOR_REAL_VIEW node_coord_z;
};

#endif  // SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
