#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BOUNDARY_CLASSIFIER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BOUNDARY_CLASSIFIER_H_

#pragma once

#include <data_type.h>
#include <model.h>

#include <cmath>

namespace model
{
/**
 * @brief Classifies boundary flags for nodes of an unstructured Cartesian mesh.
 *
 * Nodes are classified against the *global* domain bounds so that nodes on
 * interior partition boundaries (e.g. MPI subdomain edges) are not marked as
 * physical boundaries.
 *
 * Classification rules:
 *  - Not on any global face              → InteriorNode
 *  - On the z_max global face AND
 *    free_surface_on_top == true         → Surface
 *  - On any other global face            → Damping
 *
 * @tparam FloatType Floating-point type for coordinates and bounds
 * @tparam ScalarType Integer type used to cast BoundaryFlag values
 */
template <typename FloatType, typename ScalarType>
class CartesianUnstructBoundaryClassifier
{
 public:
  /**
   * @param tol Distance tolerance for boundary detection (typically
   *            min_grid_spacing * 1e-4)
   */
  CartesianUnstructBoundaryClassifier(FloatType x_min, FloatType x_max,
                                      FloatType y_min, FloatType y_max,
                                      FloatType z_min, FloatType z_max,
                                      FloatType tol, bool free_surface_on_top)
      : x_min_(x_min),
        x_max_(x_max),
        y_min_(y_min),
        y_max_(y_max),
        z_min_(z_min),
        z_max_(z_max),
        tol_(tol),
        free_surface_on_top_(free_surface_on_top)
  {
  }

  /**
   * @brief Classify every node against the global domain bounds.
   *
   * @param n_node Total number of nodes
   * @param coords_x X-coordinates of each node
   * @param coords_y Y-coordinates of each node
   * @param coords_z Z-coordinates of each node
   * @return VECTOR_INT_VIEW of size @p n_node with BoundaryFlag values
   */
  VECTOR_INT_VIEW classify(int n_node, VECTOR_REAL_VIEW coords_x,
                           VECTOR_REAL_VIEW coords_y,
                           VECTOR_REAL_VIEW coords_z) const
  {
    auto boundaries_t = allocateVector<VECTOR_INT_VIEW>(n_node, "boundaries_t");

    for (int n = 0; n < n_node; ++n)
    {
      const FloatType x = coords_x(n);
      const FloatType y = coords_y(n);
      const FloatType z = coords_z(n);

      const bool at_xmin = (fabs(x - x_min_) < tol_);
      const bool at_xmax = (fabs(x - x_max_) < tol_);
      const bool at_ymin = (fabs(y - y_min_) < tol_);
      const bool at_ymax = (fabs(y - y_max_) < tol_);
      const bool at_zmin = (fabs(z - z_min_) < tol_);
      const bool at_zmax = (fabs(z - z_max_) < tol_);

      const bool on_boundary =
          at_xmin || at_xmax || at_ymin || at_ymax || at_zmin || at_zmax;

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
