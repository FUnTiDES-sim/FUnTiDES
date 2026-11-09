#ifndef FUNTIDES_SOLVER_IMPL_SEM_LOCAL_TIME_STEP_2D_INCLUDE_SEM_SOLVER_2D_LOCAL_TIME_STEP_H_
#define FUNTIDES_SOLVER_IMPL_SEM_LOCAL_TIME_STEP_2D_INCLUDE_SEM_SOLVER_2D_LOCAL_TIME_STEP_H_

#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

#include "sem_2d_local_time_step_plan.h"
#include "sem_solver_2d.h"

namespace solver {

/**
 * @brief Experimental, source-free two-level (2:1) stabilized LF-LTS integrator.
 *
 * With lumped mass M, stiffness K and projection P onto fine nodes, the p=2
 * specialization of Grote/Michel/Sauter, arXiv:2005.13350v4, section 2.4, uses
 * B = K - a dt^2 K P M^-1 K, a = (2 delta^2 - 1)/(16 delta^2), delta = 1 + nu/4.
 * One global K application and one application on fine_support_elements suffice.
 * The recurrence is u_next = 2 u - u_prev - dt^2 M^-1 B u; the caller swaps fields.
 * Initialization and discreteEnergy() use this same symmetric modified operator.
 *
 * Rigid walls, homogeneous free surfaces and first-order absorbing boundaries are supported.
 * Diagonal damping D uses the centered derivative D (u_next-u_prev)/(2 dt).
 * Sources, sponge, adjoint and distributed runs are rejected. The plan's
 * geometric CFL estimate does NOT certify the spectrum of M^-1 B: choose a safe
 * step and validate the configuration. Reinitialize after model/geometry changes.
 *
 * Kernel launchers are public for nvcc extended-lambda compatibility. Each template
 * instantiation lives in one generated translation unit, as for SEMsolver2D.
 */
template <int ORDER, bool IS_MODEL_ON_NODES>
class SEMsolver2DLocalTimeStep : public SEMsolver2D<ORDER, IS_MODEL_ON_NODES> {
 public:
  using Base = SEMsolver2D<ORDER, IS_MODEL_ON_NODES>;
  using MeshType = typename Base::MeshType;
  using Element = typename Base::Element;
  using DataType = typename Base::DataType;
  static constexpr int kN1d = Base::kN1d;
  static constexpr int kPointsPerElement = Base::kPointsPerElement;

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size, bool surface_sponge,
                     float taper_delta) override;

  /// Select a fixed macro step and build its fine/coarse partition. nu is in [0, 0.5].
  /// nu=0 is the original, unstabilized LTS scheme; default nu=0.01 adds stabilization.
  void configureLocalTimeStep(float coarse_dt, float cfl = 0.3f, float nu = 0.01f);

  const SEM2DLocalTimeStepPlan& getLocalTimeStepPlan() const {
    requireConfigured();
    return plan_;
  }

  /// Apply K to any nodal field, with homogeneous Dirichlet rows/columns projected out.
  /// Overwrites result; input and result must not overlap. No time-step configuration needed.
  void applyStiffness(const vectorReal& field, const vectorReal& result);

  /// Apply only the element contributions in fine_support_elements. Equals the full
  /// operator on fields supported on fine_nodes; for a general field it is a partial sum.
  void applyFineSupportStiffness(const vectorReal& field, const vectorReal& result);

  /// Apply B at the configured macro step, overwriting result. Uses internal scratch.
  void applyEffectiveStiffness(const vectorReal& field, const vectorReal& result);

  /// Set the previous field at t=-dt from current u0 and optional initial velocity:
  /// u_prev = u0 - dt v0 - dt^2 M^-1 (B u0 + D v0) / 2. Empty velocity means zero velocity.
  /// Free-surface values are zeroed in both fields; the caller must not swap here.
  void initializePreviousField(DataType& data, const vectorReal& velocity = vectorReal());

  void computeForces(const float& dt, const int& timeSample, Solver::DataStruct& data) override;
  void updateSolutionForward(const float& dt, Solver::DataStruct& data) override;
  void computeOneStep(const float& dt, const int& timeSample, Solver::DataStruct& data) override {
    computeForces(dt, timeSample, data);
    updateSolutionForward(dt, data);
  }

  /// E(n-1/2) = ((u-u_prev)^T M (u-u_prev)/dt^2 + u_prev^T B u)/2.
  /// Conserved without damping; with D>=0 its loss is dt*v_center^T D v_center per step.
  /// Positivity requires the actual spectral CFL condition.
  double discreteEnergy(const DataType& data);

  /// Internal assembly launcher. elements contains validated mesh indices; all=true
  /// ignores that view. Public only because nvcc requires public kernel launchers.
  void assembleStiffness(const vectorReal& field, const vectorReal& result, const vectorInt& elements, bool all);

 private:
  void requireInitialized() const {
    if (!initialized_) throw std::logic_error("SEM2D LTS: call computeFEInit first");
  }
  void requireConfigured() const {
    requireInitialized();
    if (!configured_) throw std::logic_error("SEM2D LTS: call configureLocalTimeStep first");
  }
  void validateFields(const vectorReal& field, const vectorReal& result) const;
  void validateData(const DataType& data) const;
  void validateStep(float dt, const DataType& data) const {
    requireConfigured();
    if (dt != plan_.coarse_dt) throw std::invalid_argument("SEM2D LTS: dt must equal the configured macro step");
    validateData(data);
  }

  MeshType mesh_;
  SEM2DLocalTimeStepPlan plan_{};
  vectorInt fine_nodes_, support_elements_;
  vectorReal fine_acceleration_, correction_, energy_force_;
  float correction_factor_{0};  // a * coarse_dt^2
  bool initialized_{false};
  bool configured_{false};
};

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::computeFEInit(model::ModelApi<float, int>& mesh,
                                                                       const std::array<float, 3>& sponge_size,
                                                                       bool surface_sponge, float taper_delta) {
  initialized_ = configured_ = false;
  auto* typed = dynamic_cast<MeshType*>(&mesh);
  if (!typed || typed->getDimension() != 2 || typed->getOrder() != ORDER ||
      typed->isModelOnNodes() != IS_MODEL_ON_NODES)
    throw std::invalid_argument("SEM2D LTS: incompatible mesh, order or model location");
  for (float size : sponge_size)
    if (size != 0) throw std::invalid_argument("SEM2D LTS: sponge is not supported yet");
  if (typed->getNumberOfNodes() < 1 || typed->getNumberOfElements() < 1)
    throw std::invalid_argument("SEM2D LTS: the mesh must be nonempty");
  auto model = *typed;
  int const count = IS_MODEL_ON_NODES ? model.getNumberOfNodes() : model.getNumberOfElements();
  int invalid_material = 0;
  Kokkos::parallel_reduce(
      "SEM2D LTS validate material", Kokkos::RangePolicy<>(0, count),
      KOKKOS_LAMBDA(int const i, int& bad) {
        float const vp = IS_MODEL_ON_NODES ? model.getModelVpOnNodes(i) : model.getModelVpOnElement(i);
        float const rho = IS_MODEL_ON_NODES ? model.getModelRhoOnNodes(i) : model.getModelRhoOnElement(i);
        if (!(vp > 0) || !(rho > 0) || !Kokkos::isfinite(vp) || !Kokkos::isfinite(rho)) ++bad;
      },
      invalid_material);
  if (invalid_material) throw std::invalid_argument("SEM2D LTS: velocities and densities must be finite and positive");
  Base::computeFEInit(mesh, sponge_size, surface_sponge, taper_delta);
  auto h_mass = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), this->getMassMatrixAcoustic());
  auto h_damping = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), this->getDampingMatrix(0));
  for (size_t i = 0; i < h_mass.extent(0); ++i) {
    if (!(h_mass(i) > 0) || !std::isfinite(h_mass(i)))
      throw std::invalid_argument("SEM2D LTS: mass must be finite and positive at every node");
    if (!(h_damping(i) >= 0) || !std::isfinite(h_damping(i)))
      throw std::invalid_argument("SEM2D LTS: damping must be finite and nonnegative");
  }
  mesh_ = *typed;
  int const n = mesh_.getNumberOfNodes();
  fine_acceleration_ = vectorReal("SEM2D LTS fine acceleration", n);
  correction_ = vectorReal("SEM2D LTS correction", n);
  energy_force_ = vectorReal("SEM2D LTS energy force", n);
  initialized_ = true;
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::configureLocalTimeStep(float coarse_dt, float cfl, float nu) {
  requireInitialized();
  if (!std::isfinite(nu) || nu < 0 || nu > 0.5f)
    throw std::invalid_argument("SEM2D LTS: stabilization nu must be in [0, 0.5]");
  auto next_plan = buildSEM2DLocalTimeStepPlan(mesh_, cfl, coarse_dt);
  vectorInt nodes("SEM2D LTS fine nodes", next_plan.fine_nodes.size());
  vectorInt support("SEM2D LTS support elements", next_plan.fine_support_elements.size());
  auto h_nodes = Kokkos::create_mirror_view(nodes);
  auto h_support = Kokkos::create_mirror_view(support);
  for (size_t i = 0; i < next_plan.fine_nodes.size(); ++i) h_nodes(i) = next_plan.fine_nodes[i];
  for (size_t i = 0; i < next_plan.fine_support_elements.size(); ++i) h_support(i) = next_plan.fine_support_elements[i];
  Kokkos::deep_copy(nodes, h_nodes);
  Kokkos::deep_copy(support, h_support);
  double const delta = 1.0 + nu / 4.0;
  correction_factor_ = static_cast<float>((2 * delta * delta - 1) / (16 * delta * delta) * coarse_dt * coarse_dt);
  plan_ = std::move(next_plan);
  fine_nodes_ = nodes;
  support_elements_ = support;
  configured_ = true;
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::validateFields(const vectorReal& field,
                                                                        const vectorReal& result) const {
  requireInitialized();
  auto const n = static_cast<size_t>(mesh_.getNumberOfNodes());
  if (field.extent(0) != n || result.extent(0) != n)
    throw std::invalid_argument("SEM2D LTS: wavefields must have one entry per node");
  auto const less = std::less<const float*>();
  if (less(field.data(), result.data() + n) && less(result.data(), field.data() + n))
    throw std::invalid_argument("SEM2D LTS: input and output fields must not overlap");
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::validateData(const DataType& data) const {
  validateFields(data.getCurrentField(0), data.getPreviousField(0));
  if (data.getPrevPrevField(0).extent(0) || data.isDistributed)
    throw std::invalid_argument("SEM2D LTS: only two-buffer, single-process forward runs are supported");
  if (data.getRhsElement().extent(0) || data.getRhsTerm(0).extent(0))
    throw std::invalid_argument("SEM2D LTS: sources are not supported yet");
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::applyStiffness(const vectorReal& field,
                                                                        const vectorReal& result) {
  validateFields(field, result);
  assembleStiffness(field, result, vectorInt(), true);
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::applyFineSupportStiffness(const vectorReal& field,
                                                                                   const vectorReal& result) {
  requireConfigured();
  validateFields(field, result);
  assembleStiffness(field, result, support_elements_, false);
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::assembleStiffness(const vectorReal& field,
                                                                           const vectorReal& result,
                                                                           const vectorInt& elements, bool all) {
  Kokkos::deep_copy(result, 0.f);
  auto mesh = mesh_;
  int const count = all ? mesh.getNumberOfElements() : static_cast<int>(elements.extent(0));
  Kokkos::parallel_for(
      "SEM2D LTS stiffness", Kokkos::RangePolicy<>(0, count), KOKKOS_LAMBDA(int const item) {
        int const e = all ? item : elements(item);
        real_t u[kPointsPerElement], v[kPointsPerElement] = {0};
        int nodes[kPointsPerElement];
        for (int q = 0; q < kPointsPerElement; ++q) {
          nodes[q] = mesh.globalNodeIndex(e, q % kN1d, 0, q / kN1d);
          u[q] = mesh.isFreeSurface(nodes[q]) ? 0.f : field(nodes[q]);
        }
        real_t X[4][2];
        for (int c = 0; c < 4; ++c) {
          int const node = nodes[Element::meshIndexToLinearIndex2D(c)];
          X[c][0] = mesh.nodeCoord(node, 0);
          X[c][1] = mesh.nodeCoord(node, 2);
        }
        Element::computeStiffnessTermSumFact(X, u, v, [&](int qa, int qb) -> real_t {
          if constexpr (IS_MODEL_ON_NODES)
            return 1.f / mesh.getModelRhoOnNodes(nodes[qa + qb * kN1d]);
          else
            return 1.f / mesh.getModelRhoOnElement(e);
        });
        for (int q = 0; q < kPointsPerElement; ++q)
          if (!mesh.isFreeSurface(nodes[q])) ATOMICADD(result(nodes[q]), v[q]);
      });
  FENCE
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::applyEffectiveStiffness(const vectorReal& field,
                                                                                 const vectorReal& result) {
  requireConfigured();
  applyStiffness(field, result);
  if (plan_.fine_nodes.empty()) return;
  Kokkos::deep_copy(fine_acceleration_, 0.f);
  auto acceleration = fine_acceleration_;
  auto nodes = fine_nodes_;
  auto mass = this->getMassMatrixAcoustic();
  auto mesh = mesh_;
  Kokkos::parallel_for(
      "SEM2D LTS project fine acceleration", Kokkos::RangePolicy<>(0, nodes.extent(0)), KOKKOS_LAMBDA(int const j) {
        int const i = nodes(j);
        if (!mesh.isFreeSurface(i)) acceleration(i) = result(i) / mass(i);
      });
  FENCE
  applyFineSupportStiffness(acceleration, correction_);
  auto correction = correction_;
  float const factor = correction_factor_;
  Kokkos::parallel_for(
      "SEM2D LTS modified stiffness", Kokkos::RangePolicy<>(0, mesh.getNumberOfNodes()),
      KOKKOS_LAMBDA(int const i) { result(i) -= factor * correction(i); });
  FENCE
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::initializePreviousField(DataType& data,
                                                                                 const vectorReal& velocity) {
  requireConfigured();
  validateData(data);
  int const n = mesh_.getNumberOfNodes();
  if (velocity.extent(0) != 0 && velocity.extent(0) != static_cast<size_t>(n))
    throw std::invalid_argument("SEM2D LTS: initial velocity must be empty or have one entry per node");
  auto cur = data.getCurrentField(0);
  auto prev = data.getPreviousField(0);
  if (velocity.extent(0)) {
    validateFields(velocity, cur);
    validateFields(velocity, prev);
  }
  auto force = this->getForceVector(0);
  applyEffectiveStiffness(cur, force);
  auto mass = this->getMassMatrixAcoustic();
  auto damping = this->getDampingMatrix(0);
  auto mesh = mesh_;
  float const dt = plan_.coarse_dt;
  bool const moving = velocity.extent(0) != 0;
  Kokkos::parallel_for(
      "SEM2D LTS initial previous field", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(int const i) {
        if (mesh.isFreeSurface(i))
          cur(i) = prev(i) = 0.f;
        else
          prev(i) = cur(i) - (moving ? dt * velocity(i) : 0.f) -
                    0.5f * dt * dt * (force(i) + (moving ? damping(i) * velocity(i) : 0.f)) / mass(i);
      });
  FENCE
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::computeForces(const float& dt, const int&,
                                                                       Solver::DataStruct& data) {
  auto& acoustic = dynamic_cast<DataType&>(data);
  validateStep(dt, acoustic);
  applyEffectiveStiffness(acoustic.getCurrentField(0), this->getForceVector(0));
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::updateSolutionForward(const float& dt,
                                                                               Solver::DataStruct& data) {
  validateStep(dt, dynamic_cast<DataType&>(data));
  Base::updateSolutionForward(dt, data);
}

template <int ORDER, bool IS_MODEL_ON_NODES>
double SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>::discreteEnergy(const DataType& data) {
  requireConfigured();
  validateData(data);
  auto cur = data.getCurrentField(0);
  auto prev = data.getPreviousField(0);
  applyEffectiveStiffness(cur, energy_force_);
  auto force = energy_force_;
  auto mass = this->getMassMatrixAcoustic();
  auto mesh = mesh_;
  double const dt = plan_.coarse_dt;
  double energy = 0;
  Kokkos::parallel_reduce(
      "SEM2D LTS discrete energy", Kokkos::RangePolicy<>(0, mesh.getNumberOfNodes()),
      KOKKOS_LAMBDA(int const i, double& sum) {
        if (!mesh.isFreeSurface(i)) {
          double const v = (static_cast<double>(cur(i)) - prev(i)) / dt;
          sum += 0.5 * (mass(i) * v * v + static_cast<double>(prev(i)) * force(i));
        }
      },
      energy);
  return energy;
}

#ifndef FUNTIDES_SEM2D_LTS_INSTANTIATION
#define FUNTIDES_SEM2D_LTS_EXTERN(O)                        \
  extern template class SEMsolver2DLocalTimeStep<O, false>; \
  extern template class SEMsolver2DLocalTimeStep<O, true>;
#if MAX_SOLVER_ACOUSTIC_ORDER >= 1
FUNTIDES_SEM2D_LTS_EXTERN(1)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 2
FUNTIDES_SEM2D_LTS_EXTERN(2)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 3
FUNTIDES_SEM2D_LTS_EXTERN(3)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 4
FUNTIDES_SEM2D_LTS_EXTERN(4)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 5
FUNTIDES_SEM2D_LTS_EXTERN(5)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 6
FUNTIDES_SEM2D_LTS_EXTERN(6)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 7
FUNTIDES_SEM2D_LTS_EXTERN(7)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 8
FUNTIDES_SEM2D_LTS_EXTERN(8)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 9
FUNTIDES_SEM2D_LTS_EXTERN(9)
#endif
#undef FUNTIDES_SEM2D_LTS_EXTERN
#endif

}  // namespace solver

#endif  // FUNTIDES_SOLVER_IMPL_SEM_LOCAL_TIME_STEP_2D_INCLUDE_SEM_SOLVER_2D_LOCAL_TIME_STEP_H_
