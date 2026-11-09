#include "sem_2d_local_time_step_plan.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "gllpoints.h"

namespace solver {

SEM2DLocalTimeStepPlan buildSEM2DLocalTimeStepPlan(const model::ModelUnstruct<float, int>& mesh, float cfl,
                                                   float coarse_dt) {
  if (mesh.getDimension() != 2 || mesh.getOrder() < 1 || mesh.getOrder() > MAX_GLL_ORDER ||
      mesh.getNumberOfElements() < 1 || mesh.getNumberOfNodes() < 1)
    throw std::invalid_argument("SEM2D LTS plan: a nonempty 2-D GLL mesh is required");
  if (!(cfl > 0) || !std::isfinite(cfl) || !(coarse_dt > 0) || !std::isfinite(coarse_dt) || !(coarse_dt / 2 > 0))
    throw std::invalid_argument("SEM2D LTS plan: cfl and time steps must be finite and positive");

  int const n_elem = mesh.getNumberOfElements();
  int const n_node = mesh.getNumberOfNodes();
  int const n1d = mesh.getOrder() + 1;
  int const npe = n1d * n1d;
  if (mesh.getNumberOfPointsPerElement() != npe)
    throw std::invalid_argument("SEM2D LTS plan: inconsistent quadrilateral node count");

  arrayInt nodes("SEM2D LTS element nodes", n_elem, npe);
  arrayReal estimates("SEM2D LTS element estimates", n_elem, 3);
  Kokkos::parallel_for(
      "SEM2D LTS diagnostics", Kokkos::RangePolicy<>(0, n_elem), KOKKOS_LAMBDA(int const e) {
        bool valid = true;
        double max_vp = 0;
        for (int q = 0; q < npe; ++q) {
          int const node = mesh.globalNodeIndex(e, q % n1d, 0, q / n1d);
          nodes(e, q) = node;
          if (node < 0 || node >= n_node) {
            valid = false;
            continue;
          }
          double const vp = mesh.isModelOnNodes() ? mesh.getModelVpOnNodes(node) : mesh.getModelVpOnElement(e);
          if (!(vp > 0) || !Kokkos::isfinite(vp) || !Kokkos::isfinite(mesh.nodeCoord(node, 0)) ||
              !Kokkos::isfinite(mesh.nodeCoord(node, 2)))
            valid = false;
          if (vp > max_vp) max_vp = vp;
        }
        if (!valid) return;  // estimates is zero-initialized: rejected on the host

        double min_spacing = std::numeric_limits<double>::max();
        for (int k = 0; k < n1d; ++k)
          for (int i = 0; i < n1d; ++i) {
            int const a = nodes(e, i + k * n1d);
            for (int axis = 0; axis < 2; ++axis) {
              if ((axis == 0 ? i : k) == n1d - 1) continue;
              int const b = nodes(e, i + k * n1d + (axis == 0 ? 1 : n1d));
              double const dx = static_cast<double>(mesh.nodeCoord(b, 0)) - mesh.nodeCoord(a, 0);
              double const dz = static_cast<double>(mesh.nodeCoord(b, 2)) - mesh.nodeCoord(a, 2);
              double const spacing = Kokkos::sqrt(dx * dx + dz * dz);
              if (spacing < min_spacing) min_spacing = spacing;
            }
          }
        estimates(e, 0) = static_cast<float>(min_spacing);
        estimates(e, 1) = static_cast<float>(max_vp);
        estimates(e, 2) = static_cast<float>(cfl * min_spacing / max_vp);
      });

  auto h_estimates = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), estimates);
  auto h_nodes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nodes);
  SEM2DLocalTimeStepPlan plan;
  plan.cfl = cfl;
  plan.coarse_dt = coarse_dt;
  plan.fine_dt = coarse_dt / 2;
  plan.element_steps.reserve(n_elem);
  std::vector<bool> fine_node(n_node, false), coarse_node(n_node, false);
  double const tolerance = 1.0 + 8.0 * std::numeric_limits<float>::epsilon();
  for (int e = 0; e < n_elem; ++e) {
    float const dt = h_estimates(e, 2);
    if (!(dt > 0) || !std::isfinite(dt))
      throw std::invalid_argument("SEM2D LTS plan: invalid geometry or velocity in element " + std::to_string(e));
    if (plan.fine_dt > dt * tolerance)
      throw std::invalid_argument("SEM2D LTS plan: two levels are insufficient for element " + std::to_string(e) +
                                  "; reduce coarse_dt");
    plan.element_steps.push_back({h_estimates(e, 0), h_estimates(e, 1), dt});
    bool const fine = coarse_dt > dt * tolerance;
    (fine ? plan.fine_elements : plan.coarse_elements).push_back(e);
    for (int q = 0; q < npe; ++q) (fine ? fine_node : coarse_node)[h_nodes(e, q)] = true;
  }
  for (int n = 0; n < n_node; ++n) {
    if (fine_node[n]) {
      plan.fine_nodes.push_back(n);
      if (coarse_node[n]) plan.interface_nodes.push_back(n);
    } else {
      plan.coarse_nodes.push_back(n);
    }
  }
  for (int e = 0; e < n_elem; ++e)
    for (int q = 0; q < npe; ++q)
      if (fine_node[h_nodes(e, q)]) {
        plan.fine_support_elements.push_back(e);
        break;
      }
  return plan;
}

}  // namespace solver
