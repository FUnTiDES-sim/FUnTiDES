#ifndef FUNTIDES_SOLVER_IMPL_SEM_LOCAL_TIME_STEP_2D_INCLUDE_SEM_2D_LOCAL_TIME_STEP_PLAN_H_
#define FUNTIDES_SOLVER_IMPL_SEM_LOCAL_TIME_STEP_2D_INCLUDE_SEM_2D_LOCAL_TIME_STEP_PLAN_H_

#include <vector>

#include "model_unstruct.h"

namespace solver {

/**
 * @brief Geometry-based time-step estimate of one acoustic quadrilateral.
 *
 * estimate = cfl * min_node_spacing / max_vp. The spacing is the smallest distance
 * between adjacent GLL nodes along either reference axis. This is a diagnostic
 * estimate, not a spectral CFL bound or a stability guarantee for an LTS scheme,
 * in particular on distorted elements or across material discontinuities.
 */
struct SEM2DElementTimeStep {
  float min_node_spacing;
  float max_vp;
  float estimate;
};

/**
 * @brief Host-side partition for a future two-level SEM local-time integrator.
 *
 * Elements whose estimate is below coarse_dt use fine_dt = coarse_dt / 2.
 * Every node of a fine element is a fine degree of freedom, including nodes
 * shared with coarse elements. interface_nodes contains those shared nodes.
 * fine_support_elements contains ALL elements incident to a fine node: applying
 * the assembled operator to a field supported on fine nodes needs contributions
 * from these coarse neighbours too. Lists are sorted by mesh index.
 *
 * This structure does not advance a wavefield. The integrator must provide the
 * coupling, stability and source/boundary treatment at intermediate times.
 * Rebuild the plan after changing geometry, velocities or the time steps.
 */
struct SEM2DLocalTimeStepPlan {
  float cfl;
  float coarse_dt;
  float fine_dt;
  std::vector<SEM2DElementTimeStep> element_steps;
  std::vector<int> coarse_elements;
  std::vector<int> fine_elements;
  std::vector<int> coarse_nodes;
  std::vector<int> fine_nodes;
  std::vector<int> interface_nodes;
  std::vector<int> fine_support_elements;
};

/**
 * @brief Builds a fixed two-level (2:1) partition from a 2-D acoustic mesh.
 *
 * Reads geometry, connectivity and velocities on the device, where model updates
 * may have occurred, then returns host diagnostics and index lists. Uses the
 * maximum velocity of each element, including all its nodes for a nodal model.
 * No mass or stiffness matrix is assembled and the input model is not modified.
 *
 * @throws std::invalid_argument For a non-2-D/empty mesh, nonpositive or nonfinite
 * parameters, invalid node indices/coordinates/velocities, zero spacing, or if
 * even fine_dt exceeds an element's estimate. The latter comparison allows
 * floating-point roundoff at the level boundary.
 */
SEM2DLocalTimeStepPlan buildSEM2DLocalTimeStepPlan(const model::ModelUnstruct<float, int>& mesh, float cfl,
                                                   float coarse_dt);

}  // namespace solver

#endif  // FUNTIDES_SOLVER_IMPL_SEM_LOCAL_TIME_STEP_2D_INCLUDE_SEM_2D_LOCAL_TIME_STEP_PLAN_H_
