#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BOUNDARY_CLASSIFIER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BOUNDARY_CLASSIFIER_H_

#include <data_type.h>
#include <model.h>

#include <cmath>

namespace model {
/**
 * @brief Assigns a BoundaryFlag to each node of an unstructured Cartesian mesh.
 *
 * Nodes are tested against the global domain bounds, so nodes on interior
 * partition boundaries (for example MPI subdomain edges) are not marked as
 * physical boundaries.
 *
 * Rules:
 *  - not on any global face: InteriorNode
 *  - on the global z_max face and free_surface_on_top is true: Surface
 *  - on any other global face: Damping
 *
 * @tparam FloatType Floating-point type of coordinates and bounds.
 * @tparam ScalarType Integer type the BoundaryFlag values are cast to.
 */
template <typename FloatType, typename ScalarType>
class CartesianUnstructBoundaryClassifier {
 public:
  /**
   * @brief Builds a classifier for the given global domain.
   *
   * @param x_min Global lower bound along x.
   * @param x_max Global upper bound along x.
   * @param y_min Global lower bound along y.
   * @param y_max Global upper bound along y.
   * @param z_min Global lower bound along z.
   * @param z_max Global upper bound along z.
   * @param tol Distance below which a node is considered on a face
   *            (typically min_grid_spacing * 1e-4).
   * @param free_surface_on_top If true, nodes on the z_max face are Surface
   *                            instead of Damping.
   */
  CartesianUnstructBoundaryClassifier(FloatType x_min, FloatType x_max, FloatType y_min, FloatType y_max,
                                      FloatType z_min, FloatType z_max, FloatType tol, bool free_surface_on_top)
      : x_min_(x_min),
        x_max_(x_max),
        y_min_(y_min),
        y_max_(y_max),
        z_min_(z_min),
        z_max_(z_max),
        tol_(tol),
        free_surface_on_top_(free_surface_on_top) {}

  /**
   * @brief Classifies every node against the global domain bounds.
   *
   * @param n_node Number of nodes.
   * @param coords_x X coordinate of each node, size n_node.
   * @param coords_y Y coordinate of each node, size n_node.
   * @param coords_z Z coordinate of each node, size n_node.
   * @return Vector of size @p n_node holding one BoundaryFlag value per node.
   */
  vectorInt classify(int n_node, vectorReal coords_x, vectorReal coords_y, vectorReal coords_z) const {
    auto boundaries_t = allocateVector<vectorInt>(n_node, "boundaries_t");

    for (int n = 0; n < n_node; ++n) {
      const FloatType x = coords_x(n);
      const FloatType y = coords_y(n);
      const FloatType z = coords_z(n);

      const bool at_xmin = (fabs(x - x_min_) < tol_);
      const bool at_xmax = (fabs(x - x_max_) < tol_);
      const bool at_ymin = (fabs(y - y_min_) < tol_);
      const bool at_ymax = (fabs(y - y_max_) < tol_);
      const bool at_zmin = (fabs(z - z_min_) < tol_);
      const bool at_zmax = (fabs(z - z_max_) < tol_);

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

#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BOUNDARY_CLASSIFIER_H_
