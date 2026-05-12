#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BOUNDARY_CLASSIFIER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BOUNDARY_CLASSIFIER_H_

#include <data_type.h>
#include <model.h>

#include <cmath>

namespace model {
/**
 * @brief Classifies boundary flags for nodes of a structured Cartesian mesh.
 *
 * Node ordering follows the structured convention:
 *   n = k*(nx*ny) + j*nx + i,  with i∈[0,nx), j∈[0,ny), k∈[0,nz).
 *
 * Nodes on faces whose coordinate matches a global domain boundary (within
 * tol) receive Damping or Surface flags; all others are InteriorNode.
 *
 * Classification rules:
 *  - Not on any global face     → InteriorNode
 *  - z_max global face AND
 *    free_surface_on_top        → Surface
 *  - Any other global face      → Damping
 *
 * @tparam FloatType  Floating-point type for coordinates and bounds
 * @tparam ScalarType Integer type used to cast BoundaryFlag values
 */
template <typename FloatType, typename ScalarType>
class CartesianStructBoundaryClassifier {
 public:
  CartesianStructBoundaryClassifier(FloatType x_min, FloatType x_max, FloatType y_min, FloatType y_max, FloatType z_min,
                                    FloatType z_max, FloatType tol, bool free_surface_on_top)
      : x_min_(x_min),
        x_max_(x_max),
        y_min_(y_min),
        y_max_(y_max),
        z_min_(z_min),
        z_max_(z_max),
        tol_(tol),
        free_surface_on_top_(free_surface_on_top) {}

  /**
   * @brief Classify every node of the structured grid.
   *
   * @param n_node Total number of nodes (nx*ny*nz)
   * @param nx, ny, nz  Node counts in each dimension
   * @param ox, oy, oz  Local domain origin
   * @param lx, ly, lz  Local domain dimensions
   * @return vectorInt of size n_node with BoundaryFlag values
   */
  vectorInt classify(int n_node, int nx, int ny, int nz, FloatType ox, FloatType oy, FloatType oz, FloatType lx,
                     FloatType ly, FloatType lz) const {
    const bool x_min_is_global = fabs(ox - x_min_) < tol_;
    const bool x_max_is_global = fabs((ox + lx) - x_max_) < tol_;
    const bool y_min_is_global = fabs(oy - y_min_) < tol_;
    const bool y_max_is_global = fabs((oy + ly) - y_max_) < tol_;
    const bool z_min_is_global = fabs(oz - z_min_) < tol_;
    const bool z_max_is_global = fabs((oz + lz) - z_max_) < tol_;

    auto boundaries_t = allocateVector<vectorInt>(n_node, "boundaries_t");

    for (int n = 0; n < n_node; ++n) {
      const int i = n % nx;
      const int j = (n / nx) % ny;
      const int k = n / (nx * ny);

      const bool at_xmin = x_min_is_global && (i == 0);
      const bool at_xmax = x_max_is_global && (i == nx - 1);
      const bool at_ymin = y_min_is_global && (j == 0);
      const bool at_ymax = y_max_is_global && (j == ny - 1);
      const bool at_zmin = z_min_is_global && (k == 0);
      const bool at_zmax = z_max_is_global && (k == nz - 1);

      const bool on_boundary = at_xmin || at_xmax || at_ymin || at_ymax || at_zmin || at_zmax;

      if (!on_boundary)
        boundaries_t(n) = static_cast<ScalarType>(BoundaryFlag::InteriorNode);
      else if (free_surface_on_top_ && at_zmax)
        boundaries_t(n) = static_cast<ScalarType>(BoundaryFlag::Surface);
      else
        boundaries_t(n) = static_cast<ScalarType>(BoundaryFlag::Damping);
    }

    return boundaries_t;
  }

 private:
  FloatType x_min_, x_max_;
  FloatType y_min_, y_max_;
  FloatType z_min_, z_max_;
  FloatType tol_;
  bool free_surface_on_top_;
};

}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BOUNDARY_CLASSIFIER_H_
