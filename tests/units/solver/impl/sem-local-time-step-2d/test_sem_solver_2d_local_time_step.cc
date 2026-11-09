#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "data_type.h"
#include "gllpoints.h"
#include "sem_solver_2d_local_time_step.h"

namespace solver {
namespace test {
namespace {

// Conforming quadrilateral grid with nonuniform x widths and unit z widths.
model::ModelUnstructData<float, int> grid(int order, bool nodal, const std::vector<float>& edges, int ez = 3) {
  int const ex = edges.size() - 1, nx = ex * order + 1, nz = ez * order + 1, n1d = order + 1;
  model::ModelUnstructData<float, int> d;
  d.dim_ = 2;
  d.order_ = order;
  d.n_element_ = ex * ez;
  d.n_node_ = nx * nz;
  d.lx_ = edges.back() - edges.front();
  d.ly_ = 0;
  d.lz_ = ez;
  d.ox_ = edges.front();
  d.oy_ = d.oz_ = 0;
  d.isModelOnNodes_ = nodal;
  d.isElastic_ = false;
  d.global_node_index_ = allocateArray2D<dualArrayInt>(d.n_element_, n1d * n1d, "LTS test nodes");
  d.nodes_coords_x_ = allocateVector<dualVectorReal>(d.n_node_, "LTS test x");
  d.nodes_coords_y_ = allocateVector<dualVectorReal>(d.n_node_, "LTS test y");
  d.nodes_coords_z_ = allocateVector<dualVectorReal>(d.n_node_, "LTS test z");
  for (int k = 0; k < ez; ++k)
    for (int e = 0; e < ex; ++e)
      for (int j = 0; j < n1d; ++j)
        for (int i = 0; i < n1d; ++i) {
          int const node = e * order + i + (k * order + j) * nx;
          d.global_node_index_.view_host()(e + k * ex, i + j * n1d) = node;
          d.nodes_coords_x_.view_host()(node) =
              edges[e] + (edges[e + 1] - edges[e]) * (GLLPoints::get(order, i) + 1) / 2;
          d.nodes_coords_z_.view_host()(node) = k + (GLLPoints::get(order, j) + 1) / 2;
        }
  auto& vp = nodal ? d.model_vp_node_ : d.model_vp_element_;
  auto& rho = nodal ? d.model_rho_node_ : d.model_rho_element_;
  vp = allocateVector<dualVectorReal>(nodal ? d.n_node_ : d.n_element_, "LTS test vp");
  rho = allocateVector<dualVectorReal>(nodal ? d.n_node_ : d.n_element_, "LTS test rho");
  Kokkos::deep_copy(vp.view_host(), 1.f);
  Kokkos::deep_copy(rho.view_host(), 1.f);
  for (auto* v : {&d.nodes_coords_x_, &d.nodes_coords_y_, &d.nodes_coords_z_, &vp, &rho}) v->modify_host();
  d.global_node_index_.modify_host();
  return d;
}

std::vector<float> host(const vectorReal& v) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  return {h.data(), h.data() + h.extent(0)};
}

void fill(const vectorReal& v, const std::vector<float>& values) {
  auto h = Kokkos::create_mirror_view(v);
  for (size_t i = 0; i < values.size(); ++i) h(i) = values[i];
  Kokkos::deep_copy(v, h);
}

struct State {
  vectorReal prev, cur;
  SEMsolverDataAcoustic data;
  explicit State(int n)
      : prev("LTS test previous", n),
        cur("LTS test current", n),
        data(WavefieldAcoustic(prev, cur), RhsAcoustic(arrayReal("empty source", 0, 0), vectorInt("empty elements", 0),
                                                       arrayReal("empty weights", 0, 0))) {}
};

template <typename S>
void advance(S& solver, State& state, float dt, int count) {
  for (int n = 0; n < count; ++n) {
    solver.computeOneStep(dt, n, state.data);
    state.data.swapWavefields();
  }
}

template <int O, bool N>
struct Config {
  static constexpr int order = O;
  static constexpr bool nodal = N;
};
template <typename C>
class SEM2DLTSTest : public ::testing::Test {
 public:
  static constexpr int order = C::order;
  using LTS = SEMsolver2DLocalTimeStep<C::order, C::nodal>;
  using Global = SEMsolver2D<C::order, C::nodal>;
};
using Configs = ::testing::Types<Config<1, false>, Config<1, true>
#if MAX_SOLVER_ACOUSTIC_ORDER >= 2
                                 ,
                                 Config<2, false>, Config<2, true>
#endif
                                 >;
TYPED_TEST_SUITE(SEM2DLTSTest, Configs);

const std::vector<float> kMixedEdges{-4, -3, -2, -1, 0, 0.5f, 1, 2, 3, 4};

float macroStep(int order) { return 0.12f * (GLLPoints::get(order, 1) + 1); }

std::vector<float> smoothField(const model::ModelUnstruct<float, int>& mesh) {
  std::vector<float> values(mesh.getNumberOfNodes());
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i)
    values[i] = std::sin(0.8f * mesh.nodeCoord(i, 0)) + 0.2f * std::cos(1.3f * mesh.nodeCoord(i, 2));
  return values;
}

void absorbingEdges(model::ModelUnstructData<float, int>& d, int nx, bool free_top = false) {
  d.boundaries_t_ = allocateVector<dualVectorInt>(d.n_node_, "absorbing test boundary");
  int nz = d.n_node_ / nx;
  for (int j = 0; j < nz; ++j)
    for (int i = 0; i < nx; ++i)
      d.boundaries_t_.view_host()(i + j * nx) =
          static_cast<int>(free_top && j == nz - 1                          ? model::BoundaryFlag::Surface
                           : i == 0 || j == 0 || i == nx - 1 || j == nz - 1 ? model::BoundaryFlag::Damping
                                                                            : model::BoundaryFlag::InteriorNode);
  d.boundaries_t_.modify_host();
}

TYPED_TEST(SEM2DLTSTest, FreeTopWithAbsorbingSidesKeepsZeroPressureAndDissipationBalance) {
  auto d = grid(TestFixture::order, TypeParam::nodal, kMixedEdges);
  int nx = (kMixedEdges.size() - 1) * TestFixture::order + 1;
  absorbingEdges(d, nx, true);
  model::ModelUnstruct<float, int> mesh(d);
  typename TestFixture::LTS lts;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float dt = macroStep(TestFixture::order);
  lts.configureLocalTimeStep(dt);
  State state(mesh.getNumberOfNodes());
  std::vector<float> values(mesh.getNumberOfNodes());
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double x = (mesh.nodeCoord(i, 0) + 1) / .6, z = mesh.nodeCoord(i, 2);
    values[i] = std::exp(-.5 * x * x) *
                (std::exp(-.5 * std::pow((z - 2.5) / .3, 2)) - std::exp(-.5 * std::pow((z - 3.5) / .3, 2)));
  }
  fill(state.cur, values);
  lts.initializePreviousField(state.data);
  auto damping = host(lts.getDampingMatrix(0));
  for (int i = mesh.getNumberOfNodes() - nx; i < mesh.getNumberOfNodes(); ++i) EXPECT_FLOAT_EQ(damping[i], 0);
  ASSERT_GT(*std::max_element(damping.begin(), damping.end()), 0);
  double initial = lts.discreteEnergy(state.data), lost = 0, energy = initial;
  ASSERT_GT(initial, 0);
  for (int step = 0; step < 300; ++step) {
    auto previous = host(state.data.getPreviousField(0));
    advance(lts, state, dt, 1);
    auto current = host(state.data.getCurrentField(0));
    for (size_t i = 0; i < current.size(); ++i) {
      ASSERT_TRUE(std::isfinite(current[i]));
      if (mesh.isFreeSurface(i)) EXPECT_FLOAT_EQ(current[i], 0);
      double velocity = (static_cast<double>(current[i]) - previous[i]) / (2 * dt);
      lost += dt * damping[i] * velocity * velocity;
    }
    double next = lts.discreteEnergy(state.data);
    EXPECT_LE(next, energy + 2e-6 * initial);
    EXPECT_NEAR(next + lost, initial, 2e-5 * initial);
    energy = next;
  }
  EXPECT_LT(energy, .1 * initial);
}

TYPED_TEST(SEM2DLTSTest, AbsorbingNoFineRegionMatchesGlobalAndInitialVelocityTaylor) {
  auto d = grid(TestFixture::order, TypeParam::nodal, kMixedEdges);
  absorbingEdges(d, (kMixedEdges.size() - 1) * TestFixture::order + 1);
  model::ModelUnstruct<float, int> mesh(d);
  typename TestFixture::LTS lts;
  typename TestFixture::Global global;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  global.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float dt = macroStep(TestFixture::order) / 2;
  lts.configureLocalTimeStep(dt);
  ASSERT_TRUE(lts.getLocalTimeStepPlan().fine_nodes.empty());
  State a(mesh.getNumberOfNodes()), b(mesh.getNumberOfNodes());
  auto initial = smoothField(mesh);
  fill(a.cur, initial);
  fill(b.cur, initial);
  vectorReal velocity("initial derivative", initial.size()), stiffness("K initial", initial.size());
  auto v = initial;
  for (auto& value : v) value *= .3f;
  fill(velocity, v);
  lts.applyStiffness(a.cur, stiffness);
  auto mass = host(lts.getMassMatrixAcoustic()), damping = host(lts.getDampingMatrix(0)), force = host(stiffness);
  ASSERT_GT(*std::max_element(damping.begin(), damping.end()), 0);
  lts.initializePreviousField(a.data, velocity);
  auto previous = host(a.prev);
  for (size_t i = 0; i < previous.size(); ++i)
    EXPECT_NEAR(previous[i], initial[i] - dt * v[i] - .5f * dt * dt * (force[i] + damping[i] * v[i]) / mass[i], 1e-6);
  Kokkos::deep_copy(b.prev, a.prev);
  advance(lts, a, dt, 50);
  advance(global, b, dt, 50);
  auto actual = host(a.data.getCurrentField(0)), expected = host(b.data.getCurrentField(0));
  for (size_t i = 0; i < actual.size(); ++i) EXPECT_NEAR(actual[i], expected[i], 3e-5);
}

TYPED_TEST(SEM2DLTSTest, AbsorbingMixedLevelsSatisfyDiscreteDissipationBalance) {
  auto d = grid(TestFixture::order, TypeParam::nodal, kMixedEdges);
  absorbingEdges(d, (kMixedEdges.size() - 1) * TestFixture::order + 1);
  model::ModelUnstruct<float, int> mesh(d);
  typename TestFixture::LTS lts;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float dt = macroStep(TestFixture::order);
  lts.configureLocalTimeStep(dt);
  ASSERT_FALSE(lts.getLocalTimeStepPlan().fine_elements.empty());
  ASSERT_FALSE(lts.getLocalTimeStepPlan().coarse_elements.empty());
  State state(mesh.getNumberOfNodes());
  fill(state.cur, smoothField(mesh));
  lts.initializePreviousField(state.data);
  auto damping = host(lts.getDampingMatrix(0));
  double initial_energy = lts.discreteEnergy(state.data), energy = initial_energy, lost = 0;
  ASSERT_GT(initial_energy, 0);
  for (int step = 0; step < 500; ++step) {
    auto previous = host(state.data.getPreviousField(0));
    advance(lts, state, dt, 1);
    auto current = host(state.data.getCurrentField(0));
    for (size_t i = 0; i < current.size(); ++i) {
      ASSERT_TRUE(std::isfinite(current[i]));
      double velocity = (static_cast<double>(current[i]) - previous[i]) / (2 * dt);
      lost += dt * damping[i] * velocity * velocity;
    }
    double next_energy = lts.discreteEnergy(state.data);
    EXPECT_LE(next_energy, energy + 2e-6 * initial_energy);
    EXPECT_NEAR(next_energy + lost, initial_energy, 2e-5 * initial_energy);
    energy = next_energy;
  }
  EXPECT_LT(energy, .1 * initial_energy);
}

TYPED_TEST(SEM2DLTSTest, SpatialOperatorMatchesReferenceWithHeterogeneousDensity) {
  auto d = grid(TestFixture::order, TypeParam::nodal, kMixedEdges);
  auto& rho = TypeParam::nodal ? d.model_rho_node_ : d.model_rho_element_;
  for (size_t i = 0; i < rho.extent(0); ++i) rho.view_host()(i) = 1.f + 0.03f * (i % 7);
  model::ModelUnstruct<float, int> mesh(d);
  typename TestFixture::LTS lts;
  typename TestFixture::Global global;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  global.computeFEInit(mesh, {0, 0, 0}, false, 0);
  State state(mesh.getNumberOfNodes());
  fill(state.cur, smoothField(mesh));
  vectorReal result("spatial output", mesh.getNumberOfNodes());
  lts.applyStiffness(state.cur, result);
  global.computeForces(0.01f, 0, state.data);
  auto expected = host(global.getForceVector(0)), actual = host(result);
  for (size_t i = 0; i < actual.size(); ++i) EXPECT_NEAR(actual[i], expected[i], 3e-6);
}

TYPED_TEST(SEM2DLTSTest, FineSupportMatchesFullOperatorOnFineFields) {
  model::ModelUnstruct<float, int> mesh(grid(TestFixture::order, TypeParam::nodal, kMixedEdges));
  typename TestFixture::LTS lts;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  lts.configureLocalTimeStep(macroStep(TestFixture::order));
  auto const& plan = lts.getLocalTimeStepPlan();
  ASSERT_FALSE(plan.fine_elements.empty());
  ASSERT_GT(plan.fine_support_elements.size(), plan.fine_elements.size());
  ASSERT_LT(plan.fine_support_elements.size(), static_cast<size_t>(mesh.getNumberOfElements()));
  vectorReal field("fine field", mesh.getNumberOfNodes()), full("full K", mesh.getNumberOfNodes()),
      partial("local K", mesh.getNumberOfNodes());
  auto values = smoothField(mesh);
  std::vector<float> restricted(values.size(), 0.f);
  for (int i : plan.fine_nodes) restricted[i] = values[i];
  fill(field, restricted);
  lts.applyStiffness(field, full);
  lts.applyFineSupportStiffness(field, partial);
  auto a = host(full), b = host(partial);
  for (size_t i = 0; i < a.size(); ++i) EXPECT_NEAR(a[i], b[i], 2e-6);
  // A general field must still see a genuinely restricted operator.
  fill(field, values);
  lts.applyStiffness(field, full);
  lts.applyFineSupportStiffness(field, partial);
  a = host(full);
  b = host(partial);
  double difference = 0;
  for (size_t i = 0; i < a.size(); ++i) difference += std::abs(a[i] - b[i]);
  EXPECT_GT(difference, 0.1);
}

TYPED_TEST(SEM2DLTSTest, NoFineRegionMatchesGlobalVerlet) {
  model::ModelUnstruct<float, int> mesh(grid(TestFixture::order, TypeParam::nodal, kMixedEdges));
  typename TestFixture::LTS lts;
  typename TestFixture::Global global;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  global.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float const dt = macroStep(TestFixture::order) / 2;
  lts.configureLocalTimeStep(dt);
  ASSERT_TRUE(lts.getLocalTimeStepPlan().fine_elements.empty());
  State a(mesh.getNumberOfNodes()), b(mesh.getNumberOfNodes());
  auto values = smoothField(mesh);
  fill(a.cur, values);
  fill(b.cur, values);
  lts.initializePreviousField(a.data);
  Kokkos::deep_copy(b.prev, a.prev);
  advance(lts, a, dt, 50);
  advance(global, b, dt, 50);
  auto actual = host(a.data.getCurrentField(0)), expected = host(b.data.getCurrentField(0));
  for (size_t i = 0; i < actual.size(); ++i) EXPECT_NEAR(actual[i], expected[i], 3e-5);
}

TYPED_TEST(SEM2DLTSTest, MatchesIndependentTwoStageChebyshevRecurrence) {
  model::ModelUnstruct<float, int> mesh(grid(TestFixture::order, TypeParam::nodal, kMixedEdges));
  typename TestFixture::LTS lts;
  typename TestFixture::Global global;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  global.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float const dt = macroStep(TestFixture::order), nu = 0.3f;
  lts.configureLocalTimeStep(dt, 0.3f, nu);
  State state(mesh.getNumberOfNodes()), reference(mesh.getNumberOfNodes());
  auto u = smoothField(mesh), previous = u;
  for (auto& x : previous) x *= 0.9f;
  fill(state.cur, u);
  fill(state.prev, previous);
  auto mass = host(global.getMassMatrixAcoustic());
  auto applyA = [&](const std::vector<float>& v) {
    fill(reference.cur, v);
    global.computeForces(dt, 0, reference.data);
    auto result = host(global.getForceVector(0));
    for (size_t i = 0; i < result.size(); ++i) result[i] /= mass[i];
    return result;
  };
  auto coarse = u;
  for (int i : lts.getLocalTimeStepPlan().fine_nodes) coarse[i] = 0;
  auto w = applyA(coarse), au = applyA(u), z1 = u;
  double const delta = 1 + nu / 4.0, t2 = 2 * delta * delta - 1, omega = 8 * delta / t2;
  for (size_t i = 0; i < u.size(); ++i) z1[i] -= dt * dt / (omega * delta) * au[i];
  std::vector<float> fine(u.size(), 0);
  for (int i : lts.getLocalTimeStepPlan().fine_nodes) fine[i] = z1[i];
  auto afine = applyA(fine);
  lts.computeOneStep(dt, 0, state.data);
  auto actual = host(state.data.getPreviousField(0));
  for (size_t i = 0; i < u.size(); ++i) {
    double const z2 = (1 + 1 / t2) * z1[i] - u[i] / t2 - 2 * dt * dt / omega * delta / t2 * (w[i] + afine[i]);
    EXPECT_NEAR(actual[i], 2 * z2 - previous[i], 2e-6);
  }
}

TYPED_TEST(SEM2DLTSTest, ModifiedEnergyIsConservedAcrossTimeInterface) {
  model::ModelUnstruct<float, int> mesh(grid(TestFixture::order, TypeParam::nodal, kMixedEdges));
  typename TestFixture::LTS lts;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float const dt = macroStep(TestFixture::order);
  lts.configureLocalTimeStep(dt);
  ASSERT_FALSE(lts.getLocalTimeStepPlan().interface_nodes.empty());
  State state(mesh.getNumberOfNodes());
  fill(state.cur, smoothField(mesh));
  lts.initializePreviousField(state.data);
  double const initial = lts.discreteEnergy(state.data);
  ASSERT_GT(initial, 0);
  double max_drift = 0;
  for (int n = 0; n < 1500; ++n) {
    advance(lts, state, dt, 1);
    if (n % 50 == 0) max_drift = std::max(max_drift, std::abs(lts.discreteEnergy(state.data) / initial - 1));
  }
  std::cout << "[LTS energy] order=" << TestFixture::order << " nodal=" << TypeParam::nodal
            << " max relative drift=" << max_drift << '\n';
  EXPECT_LT(max_drift, 2e-4);
}

TYPED_TEST(SEM2DLTSTest, RejectsInvalidMaterialBeforeAssembly) {
  for (float invalid : {0.f, -1.f, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()}) {
    for (bool change_speed : {false, true}) {
      auto d = grid(TestFixture::order, TypeParam::nodal, kMixedEdges);
      auto& vp = TypeParam::nodal ? d.model_vp_node_ : d.model_vp_element_;
      auto& rho = TypeParam::nodal ? d.model_rho_node_ : d.model_rho_element_;
      (change_speed ? vp : rho).view_host()(0) = invalid;
      model::ModelUnstruct<float, int> mesh(d);
      typename TestFixture::LTS lts;
      EXPECT_THROW(lts.computeFEInit(mesh, {0, 0, 0}, false, 0), std::invalid_argument);
    }
  }
}

TEST(SEM2DLTS, RejectsUnsupportedRunsAndInvalidState) {
  model::ModelUnstruct<float, int> mesh(grid(1, false, kMixedEdges));
  SEMsolver2DLocalTimeStep<1, false> lts;
  State state(mesh.getNumberOfNodes());
  EXPECT_THROW(lts.configureLocalTimeStep(0.1f), std::logic_error);
  EXPECT_THROW(lts.computeFEInit(mesh, {1, 0, 0}, false, 0), std::invalid_argument);
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  EXPECT_THROW(lts.computeOneStep(0.1f, 0, state.data), std::logic_error);
  EXPECT_THROW(lts.configureLocalTimeStep(0.24f, 0.3f, -0.01f), std::invalid_argument);
  EXPECT_THROW(lts.configureLocalTimeStep(0.24f, 0.3f, 0.51f), std::invalid_argument);
  EXPECT_THROW(lts.configureLocalTimeStep(0.24f, 0.3f, std::numeric_limits<float>::quiet_NaN()), std::invalid_argument);
  lts.configureLocalTimeStep(0.24f);
  EXPECT_THROW(lts.computeOneStep(0.12f, 0, state.data), std::invalid_argument);
  state.data.isDistributed = true;
  EXPECT_THROW(lts.computeOneStep(0.24f, 0, state.data), std::invalid_argument);
  state.data.isDistributed = false;
  SEMsolverDataAcoustic sourced(
      WavefieldAcoustic(state.prev, state.cur),
      RhsAcoustic(arrayReal("source", 1, 1), vectorInt("element", 1), arrayReal("weights", 1, 4)));
  EXPECT_THROW(lts.computeOneStep(0.24f, 0, sourced), std::invalid_argument);
  SEMsolverDataAcoustic adjoint(WavefieldAcoustic(vectorReal("third", mesh.getNumberOfNodes()), state.prev, state.cur),
                                RhsAcoustic(arrayReal("none", 0, 0), vectorInt("none", 0), arrayReal("none", 0, 0)));
  EXPECT_THROW(lts.computeOneStep(0.24f, 0, adjoint), std::invalid_argument);
  EXPECT_THROW(lts.applyStiffness(state.cur, state.cur), std::invalid_argument);
  vectorReal wrong("wrong size", 1);
  EXPECT_THROW(lts.applyStiffness(state.cur, wrong), std::invalid_argument);
  EXPECT_THROW(lts.initializePreviousField(state.data, wrong), std::invalid_argument);
  EXPECT_THROW(lts.initializePreviousField(state.data, state.cur), std::invalid_argument);
  // Reinitialization invalidates the previous temporal configuration.
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  EXPECT_THROW(lts.computeOneStep(0.24f, 0, state.data), std::logic_error);
}

TEST(SEM2DLTS, AcceptsAbsorbingBoundaryAndRejectsModelLocationMismatch) {
  auto d = grid(1, false, kMixedEdges);
  d.boundaries_t_ = allocateVector<dualVectorInt>(d.n_node_, "damping boundary");
  d.boundaries_t_.view_host()(0) = static_cast<int>(model::BoundaryFlag::Damping);
  d.boundaries_t_.modify_host();
  model::ModelUnstruct<float, int> damped(d);
  SEMsolver2DLocalTimeStep<1, false> lts;
  EXPECT_NO_THROW(lts.computeFEInit(damped, {0, 0, 0}, false, 0));
  model::ModelUnstruct<float, int> nodal(grid(1, true, kMixedEdges));
  EXPECT_THROW(lts.computeFEInit(nodal, {0, 0, 0}, false, 0), std::invalid_argument);
}

TEST(SEM2DLTS, FreeSurfaceProjectionRemainsSymmetricAndZero) {
  auto d = grid(1, false, kMixedEdges);
  d.boundaries_t_ = allocateVector<dualVectorInt>(d.n_node_, "free surface");
  int const nx = kMixedEdges.size();
  for (int i = d.n_node_ - nx; i < d.n_node_; ++i)
    d.boundaries_t_.view_host()(i) = static_cast<int>(model::BoundaryFlag::Surface);
  d.boundaries_t_.modify_host();
  model::ModelUnstruct<float, int> mesh(d);
  SEMsolver2DLocalTimeStep<1, false> lts;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  lts.configureLocalTimeStep(0.24f);
  State state(d.n_node_);
  fill(state.cur, smoothField(mesh));
  lts.initializePreviousField(state.data);
  double const energy = lts.discreteEnergy(state.data);
  advance(lts, state, 0.24f, 300);
  EXPECT_NEAR(lts.discreteEnergy(state.data) / energy, 1.0, 1e-4);
  auto cur = host(state.data.getCurrentField(0)), prev = host(state.data.getPreviousField(0));
  for (int i = d.n_node_ - nx; i < d.n_node_; ++i) {
    EXPECT_EQ(cur[i], 0);
    EXPECT_EQ(prev[i], 0);
  }
}

#if MAX_SOLVER_ACOUSTIC_ORDER >= 2
TEST(SEM2DLTS, CornerNeighboursContributeToFineOperator) {
  auto d = grid(2, false, {-2, -1, 0, 1, 2, 3}, 5);
  d.model_vp_element_.view_host()(12) = 2.f;
  model::ModelUnstruct<float, int> mesh(d);
  SEMsolver2DLocalTimeStep<2, false> lts;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  lts.configureLocalTimeStep(0.12f);
  auto const& plan = lts.getLocalTimeStepPlan();
  EXPECT_EQ(plan.fine_elements, (std::vector<int>{12}));
  EXPECT_EQ(plan.fine_support_elements, (std::vector<int>{6, 7, 8, 11, 12, 13, 16, 17, 18}));
  EXPECT_EQ(plan.interface_nodes.size(), 8u);
  vectorReal field("fine corner field", d.n_node_), full("full corner result", d.n_node_),
      partial("local corner result", d.n_node_);
  std::vector<float> values(d.n_node_, 0.f);
  for (int i : plan.fine_nodes) values[i] = 1.f + 0.1f * (i % 7);
  fill(field, values);
  lts.applyStiffness(field, full);
  lts.applyFineSupportStiffness(field, partial);
  auto a = host(full), b = host(partial);
  for (size_t i = 0; i < a.size(); ++i) EXPECT_NEAR(a[i], b[i], 2e-6);
}

// A right-going, z-independent pulse crosses both time interfaces. The reference
// uses the same spatial discretization and a small global step, avoiding spatial
// error contamination in the temporal convergence measurement.
TEST(SEM2DLTS, PulseCrossesInterfacesWithSecondOrderTemporalConvergence) {
  std::vector<float> edges{-8, -7, -6, -5, -4, -3, -2, -1, 0, 0.5f, 1, 1.5f, 2, 3, 4, 5, 6, 7, 8};
  model::ModelUnstruct<float, int> mesh(grid(2, false, edges, 1));
  int const nn = mesh.getNumberOfNodes();
  float const coarse_dt = 0.12f;
  int const coarse_steps = 50;
  std::vector<float> initial(nn), velocity(nn);
  for (int i = 0; i < nn; ++i) {
    float const x = mesh.nodeCoord(i, 0) + 3;
    initial[i] = std::exp(-x * x / 2);
    velocity[i] = x * initial[i];  // -d/dx exp(-(x+3)^2/2)
  }
  SEMsolver2D<2, false> reference;
  reference.computeFEInit(mesh, {0, 0, 0}, false, 0);
  auto mass = host(reference.getMassMatrixAcoustic());
  State ref(nn);
  fill(ref.cur, initial);
  float const ref_dt = coarse_dt / 16;
  reference.computeForces(ref_dt, 0, ref.data);
  auto force = host(reference.getForceVector(0)), previous = initial;
  for (int i = 0; i < nn; ++i) previous[i] -= ref_dt * velocity[i] + 0.5f * ref_dt * ref_dt * force[i] / mass[i];
  fill(ref.prev, previous);
  advance(reference, ref, ref_dt, coarse_steps * 16);
  auto expected = host(ref.data.getCurrentField(0));
  std::vector<double> errors;
  vectorReal v("initial pulse velocity", nn);
  fill(v, velocity);
  for (int refinement : {1, 2, 4}) {
    SEMsolver2DLocalTimeStep<2, false> lts;
    lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
    float const dt = coarse_dt / refinement;
    lts.configureLocalTimeStep(dt, 0.3f / refinement);
    ASSERT_EQ(lts.getLocalTimeStepPlan().fine_elements.size(), 4u);
    State state(nn);
    fill(state.cur, initial);
    lts.initializePreviousField(state.data, v);
    advance(lts, state, dt, coarse_steps * refinement);
    auto actual = host(state.data.getCurrentField(0));
    double err = 0, norm = 0;
    for (int i = 0; i < nn; ++i) {
      err += mass[i] * std::pow(actual[i] - expected[i], 2);
      norm += mass[i] * expected[i] * expected[i];
    }
    errors.push_back(std::sqrt(err / norm));
  }
  // Float32 floor of the time integration: global Verlet at the finest fine step
  // (coarse_dt / 8). About 4e-4 on GPU (FMA rounding), well below on CPU.
  double floor_error = 0;
  {
    SEMsolver2D<2, false> global;
    global.computeFEInit(mesh, {0, 0, 0}, false, 0);
    float const dt = coarse_dt / 8;
    State state(nn);
    fill(state.cur, initial);
    global.computeForces(dt, 0, state.data);
    auto f0 = host(global.getForceVector(0));
    std::vector<float> start = initial;
    for (int i = 0; i < nn; ++i) start[i] -= dt * velocity[i] + 0.5f * dt * dt * f0[i] / mass[i];
    fill(state.prev, start);
    advance(global, state, dt, coarse_steps * 8);
    auto actual = host(state.data.getCurrentField(0));
    double err = 0, norm = 0;
    for (int i = 0; i < nn; ++i) {
      err += mass[i] * std::pow(actual[i] - expected[i], 2);
      norm += mass[i] * expected[i] * expected[i];
    }
    floor_error = std::sqrt(err / norm);
  }
  std::cout << "[LTS pulse] relative errors: " << errors[0] << ", " << errors[1] << ", " << errors[2]
            << " (global Verlet at the finest fine step: " << floor_error << ")\n";
  EXPECT_LT(errors[0], 0.03);
  EXPECT_GT(errors[0] / errors[1], 3.0);
  // Second order down to the float32 floor of the time integration
  if (errors[2] > 1.5 * floor_error) EXPECT_GT(errors[1] / errors[2], 3.0);
}

TEST(SEM2DLTS, RemainsBoundedWhenGlobalMacroStepIsUnstable) {
  std::vector<float> edges{-5, -4, -3, -2, -1, 0, 0.5f, 1, 1.5f, 2, 3, 4, 5};
  model::ModelUnstruct<float, int> mesh(grid(2, false, edges, 6));
  SEMsolver2DLocalTimeStep<2, false> lts;
  SEMsolver2D<2, false> global;
  lts.computeFEInit(mesh, {0, 0, 0}, false, 0);
  global.computeFEInit(mesh, {0, 0, 0}, false, 0);
  float const dt = 0.25f;
  lts.configureLocalTimeStep(dt, 0.8f);
  ASSERT_FALSE(lts.getLocalTimeStepPlan().coarse_elements.empty());
  ASSERT_FALSE(lts.getLocalTimeStepPlan().fine_elements.empty());
  State a(mesh.getNumberOfNodes()), b(mesh.getNumberOfNodes());
  auto values = smoothField(mesh);
  for (size_t i = 0; i < values.size(); ++i) values[i] += 0.1f * std::sin(static_cast<float>(i * 17));
  fill(a.cur, values);
  fill(b.cur, values);
  lts.initializePreviousField(a.data);
  Kokkos::deep_copy(b.prev, a.prev);
  double const energy = lts.discreteEnergy(a.data);
  ASSERT_GT(energy, 0);
  advance(lts, a, dt, 500);
  EXPECT_NEAR(lts.discreteEnergy(a.data) / energy, 1.0, 2e-4);
  advance(global, b, dt, 50);
  double max_global = 0;
  for (float x : host(b.data.getCurrentField(0))) max_global = std::max(max_global, std::abs(static_cast<double>(x)));
  EXPECT_GT(max_global, 1e5);
}
#endif

}  // namespace
}  // namespace test
}  // namespace solver
