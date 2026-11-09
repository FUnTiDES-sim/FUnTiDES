#ifndef FUNTIDES_SEM2D_LTS_VALIDATION_THREE_LAYER_CASE_H_
#define FUNTIDES_SEM2D_LTS_VALIDATION_THREE_LAYER_CASE_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "data_type.h"
#include "gllpoints.h"
#include "sem_solver_2d_local_time_step.h"

namespace solver::validation {

struct ThreeLayerSnapshot {
  float time;
  std::vector<float> local, reference;
};

struct ThreeLayerHistory {
  float time;
  double energy, relative_error;
};

struct ThreeLayerResult {
  static constexpr int ex = 60, ez = 30, order = 2, nx = ex * order + 1, nz = ez * order + 1;
  static constexpr float h = 100.f, coarse_dt = 0.012f, reference_dt = coarse_dt / 8, cfl = 0.4f;
  std::vector<float> x, z, vp;
  SEM2DLocalTimeStepPlan plan;
  std::vector<ThreeLayerSnapshot> snapshots;
  std::vector<ThreeLayerHistory> history;
};

inline model::ModelUnstructData<float, int> threeLayerMesh() {
  using C = ThreeLayerResult;
  model::ModelUnstructData<float, int> d;
  d.dim_ = 2;
  d.order_ = C::order;
  d.n_element_ = C::ex * C::ez;
  d.n_node_ = C::nx * C::nz;
  d.lx_ = C::ex * C::h;
  d.ly_ = 0;
  d.lz_ = C::ez * C::h;
  d.ox_ = d.oy_ = d.oz_ = 0;
  d.isModelOnNodes_ = false;
  d.isElastic_ = false;
  d.global_node_index_ = allocateArray2D<dualArrayInt>(d.n_element_, 9, "three layer nodes");
  d.nodes_coords_x_ = allocateVector<dualVectorReal>(d.n_node_, "three layer x");
  d.nodes_coords_y_ = allocateVector<dualVectorReal>(d.n_node_, "three layer y");
  d.nodes_coords_z_ = allocateVector<dualVectorReal>(d.n_node_, "three layer z");
  d.model_vp_element_ = allocateVector<dualVectorReal>(d.n_element_, "three layer vp");
  d.model_rho_element_ = allocateVector<dualVectorReal>(d.n_element_, "three layer rho");
  for (int k = 0; k < C::ez; ++k)
    for (int e = 0; e < C::ex; ++e) {
      int const element = e + k * C::ex;
      // z points upward: fastest layer at the bottom, slowest at the top.
      d.model_vp_element_.view_host()(element) = k < 10 ? 3000.f : k < 20 ? 2250.f : 1500.f;
      d.model_rho_element_.view_host()(element) = 1.f;
      for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
          int const node = e * 2 + i + (k * 2 + j) * C::nx;
          d.global_node_index_.view_host()(element, i + j * 3) = node;
          d.nodes_coords_x_.view_host()(node) = C::h * (e + (GLLPoints::get(2, i) + 1) / 2);
          d.nodes_coords_z_.view_host()(node) = C::h * (k + (GLLPoints::get(2, j) + 1) / 2);
        }
    }
  for (auto* v :
       {&d.nodes_coords_x_, &d.nodes_coords_y_, &d.nodes_coords_z_, &d.model_vp_element_, &d.model_rho_element_})
    v->modify_host();
  d.global_node_index_.modify_host();
  return d;
}

inline std::vector<float> layerHost(const vectorReal& field) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field);
  return {h.data(), h.data() + h.extent(0)};
}

// The test and figure exporter run this exact same numerical experiment.
inline ThreeLayerResult runThreeLayerValidation() {
  using C = ThreeLayerResult;
  auto d = threeLayerMesh();
  model::ModelUnstruct<float, int> mesh(d);
  SEMsolver2DLocalTimeStep<2, false> lts;
  SEMsolver2D<2, false> reference;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  lts.configureLocalTimeStep(C::coarse_dt, C::cfl);
  reference.computeFEInit(mesh, {0, 0, 0}, false, 0);

  ThreeLayerResult result;
  result.plan = lts.getLocalTimeStepPlan();
  result.x.resize(d.n_node_);
  result.z.resize(d.n_node_);
  for (int i = 0; i < d.n_node_; ++i) {
    result.x[i] = mesh.nodeCoord(i, 0);
    result.z[i] = mesh.nodeCoord(i, 2);
  }
  result.vp.resize(d.n_element_);
  for (int e = 0; e < d.n_element_; ++e) result.vp[e] = mesh.getModelVpOnElement(e);

  vectorReal cur("layers local current", d.n_node_), prev("layers local previous", d.n_node_);
  vectorReal ref_cur("layers reference current", d.n_node_), ref_prev("layers reference previous", d.n_node_);
  auto initial = Kokkos::create_mirror_view(cur);
  for (int i = 0; i < d.n_node_; ++i) {
    double const dx = (result.x[i] - 3000) / 160., dz = (result.z[i] - 2500) / 160.;
    initial(i) = static_cast<float>(std::exp(-0.5 * (dx * dx + dz * dz)));
  }
  Kokkos::deep_copy(cur, initial);
  Kokkos::deep_copy(ref_cur, initial);
  RhsAcoustic empty(arrayReal("layers no source", 0, 0), vectorInt("layers no element", 0),
                    arrayReal("layers no weights", 0, 0));
  SEMsolverDataAcoustic data(WavefieldAcoustic(prev, cur), empty);
  SEMsolverDataAcoustic ref_data(WavefieldAcoustic(ref_prev, ref_cur), empty);
  lts.initializePreviousField(data);  // zero initial time derivative
  reference.computeForces(C::reference_dt, 0, ref_data);
  auto force = layerHost(reference.getForceVector(0));
  auto mass = layerHost(reference.getMassMatrixAcoustic());
  auto ref_start = Kokkos::create_mirror_view(ref_prev);
  for (int i = 0; i < d.n_node_; ++i)
    ref_start(i) = initial(i) - 0.5f * C::reference_dt * C::reference_dt * force[i] / mass[i];
  Kokkos::deep_copy(ref_prev, ref_start);

  auto record = [&](int step, bool snapshot) {
    auto a = layerHost(data.getCurrentField(0)), b = layerHost(ref_data.getCurrentField(0));
    double error = 0, norm = 0;
    for (int i = 0; i < d.n_node_; ++i) {
      error += static_cast<double>(mass[i]) * (a[i] - b[i]) * (a[i] - b[i]);
      norm += static_cast<double>(mass[i]) * b[i] * b[i];
    }
    result.history.push_back({step * C::coarse_dt, lts.discreteEnergy(data), std::sqrt(error / norm)});
    if (snapshot) result.snapshots.push_back({step * C::coarse_dt, std::move(a), std::move(b)});
  };
  record(0, true);
  for (int step = 1; step <= 150; ++step) {
    lts.computeOneStep(C::coarse_dt, step - 1, data);
    data.swapWavefields();
    for (int j = 0; j < 8; ++j) {
      reference.computeOneStep(C::reference_dt, (step - 1) * 8 + j, ref_data);
      ref_data.swapWavefields();
    }
    record(step, step == 30 || step == 60 || step == 90 || step == 150);
  }
  return result;
}

}  // namespace solver::validation

#endif
