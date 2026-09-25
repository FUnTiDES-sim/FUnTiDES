#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BOUNDARY_CLASSIFIER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BOUNDARY_CLASSIFIER_H_

#include <data_type.h>
#include <model.h>

#include <cmath>

namespace model {
/**
 * @brief Assigns a BoundaryFlag to every node of a structured Cartesian subdomain.
 *
 * A node is on a global boundary when it lies on a face of the local subdomain
 * whose coordinate matches a face of the global domain within a tolerance.
 * Such nodes get Surface (top face with a free surface) or Damping; all other
 * nodes get InteriorNode. Faces shared with a neighbouring subdomain are
 * never boundaries.
 *
 * Node numbering is n = k*(nx*ny) + j*nx + i, with i in [0,nx), j in [0,ny),
 * k in [0,nz).
 *
 * @tparam FloatType  Floating-point type of coordinates and bounds.
 * @tparam ScalarType Integer type the BoundaryFlag values are cast to.
 */
template <typename FloatType, typename ScalarType>
class CartesianStructBoundaryClassifier {
 public:
  /**
   * @brief Stores the global domain bounds and the classification options.
   *
   * @param x_min,x_max,y_min,y_max,z_min,z_max Global domain bounds.
   * @param tol Tolerance used to compare subdomain faces with the global bounds.
   * @param free_surface_on_top If true, nodes on the global z_max face are Surface instead of Damping.
   * @todo VERIFY: are the bounds and tol in the same length unit as the origin and size passed to classify()?
   */
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
   * @brief Classifies every node of the local structured grid.
   *
   * @param n_node Total number of nodes, must equal nx*ny*nz.
   * @param nx,ny,nz Node counts along each axis.
   * @param ox,oy,oz Origin of the local subdomain.
   * @param lx,ly,lz Extent of the local subdomain along each axis.
   * @return Newly allocated vector of size n_node, indexed by node number, holding BoundaryFlag values.
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
