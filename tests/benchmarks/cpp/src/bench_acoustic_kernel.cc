// =============================================================================
// Suggested run:
//   ./bench_acoustic_kernels \
//       --benchmark_min_time=2s \
//       --benchmark_repetitions=10 \
//       --benchmark_report_aggregates_only=true
// =============================================================================

#include <algorithm>
#include <cmath>

#include "Integrals.h"
#include "Qk_Hexahedron_Tensorial.h"
#include "bench_main.h"
#include "bench_solver_fe_cartesian_unstruct_common.h"
#include "model_unstruct.h"
#include "rhs_acoustic.h"
#include "sem_enums.h"
#include "sem_solver.h"
#include "sem_solver_impl.h"
#include "wavefield_acoustic.h"

template <int ORDER>
using Qk_Hexahedron_Tensorial_Selector = Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<ORDER>;

namespace model {
namespace bench {

template <int Order>
using AcousticKernelConfig = BuilderConfig<Order, false>;

template <typename T>
class AcousticKernelFixture : public SolverUnstructFixture<T> {};

struct BenchmarkArrays {
  arrayReal rhsTerm;
  vectorInt rhsElement;
  arrayReal rhsWeights;
  vectorReal pnGlobalPrev;
  vectorReal pnGlobalCurr;
  arrayReal rhsLocation;

  BenchmarkArrays(int n_rhs, int n_time_steps, int n_dof, int nb_points_per_element) {
    rhsTerm = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTerm");
    rhsElement = allocateVector<vectorInt>(n_rhs, "rhsElement");
    rhsWeights = allocateArray2D<arrayReal>(n_rhs, nb_points_per_element, "rhsWeights");
    pnGlobalPrev = allocateVector<vectorReal>(n_dof, "pnGlobalPrev");
    pnGlobalCurr = allocateVector<vectorReal>(n_dof, "pnGlobalCurr");
    rhsLocation = allocateArray2D<arrayReal>(1, 3, "rhsLocation");
    FENCE
  }
};

namespace {

// Time a single device call: fence, start, call, fence, stop. Returns seconds.
template <typename Fn>
inline double timeDevice(Fn&& fn) {
  Kokkos::fence();
  Kokkos::Timer timer;
  fn();
  Kokkos::fence();
  return timer.seconds();
}

// DOF-updates per second + the problem size as a visible counter.
inline void setThroughput(benchmark::State& state, int64_t n_dof) {
  state.SetItemsProcessed(state.iterations() * n_dof);  // -> items_per_second == DOF/s
  state.counters["nDOF"] = static_cast<double>(n_dof);
}

// Max relative difference between two solvers' assembled force vectors.
template <typename SolverRef, typename SolverTest>
inline double maxRelDiffForce(SolverRef& ref, SolverTest& test, int n_fields) {
  double dmax = 0.0, rmax = 0.0;
  for (int f = 0; f < n_fields; ++f) {
    auto hr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, ref.getForceVector(f));
    auto ht = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, test.getForceVector(f));
    const int n = static_cast<int>(hr.extent(0));
    for (int i = 0; i < n; ++i) {
      dmax = std::max(dmax, std::abs(static_cast<double>(hr(i)) - static_cast<double>(ht(i))));
      rmax = std::max(rmax, std::abs(static_cast<double>(hr(i))));
    }
  }
  return rmax > 0.0 ? dmax / rmax : dmax;
}

}  // namespace

// ==============================================================================
// BENCHMARK 1 : Flat (makutu)
// ==============================================================================
BENCHMARK_TEMPLATE_METHOD_F(AcousticKernelFixture, ElementContrib_Flat)(benchmark::State& state) {
  auto model = this->createModel();
  auto solver =
      solver_factory::createSolver(methodType::kSem, this->implem_, meshType::kUnstruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kAcoustic, this->order);
  solver->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);

  BenchmarkArrays arrays(this->n_rhs, this->n_time_steps, this->n_dof, model->getNumberOfPointsPerElement());
  auto wavefield = WavefieldAcoustic(arrays.pnGlobalPrev, arrays.pnGlobalCurr);
  auto rhs = RhsAcoustic(arrays.rhsTerm, arrays.rhsElement, arrays.rhsWeights);
  SEMsolverDataAcoustic data(wavefield, rhs);

  using Integral = typename IntegralTypeSelector<this->order, IntegralType::MAKUTU>::type;
  using Model = ModelUnstruct<float, int>;
  using SolverOnNodes = solver::fe::SEMsolver<this->order, Integral, Model, true, solver::fe::physicType::kAcoustic>;
  using SolverOnElems = solver::fe::SEMsolver<this->order, Integral, Model, false, solver::fe::physicType::kAcoustic>;

  if (this->isModelOnNodes_) {
    auto* s = dynamic_cast<SolverOnNodes*>(solver.get());
    s->computeElementContributions_Acoustic_Flat(data);  // warmup
    FENCE
    for (auto _ : state)
      state.SetIterationTime(timeDevice([&] { s->computeElementContributions_Acoustic_Flat(data); }));
  } else {
    auto* s = dynamic_cast<SolverOnElems*>(solver.get());
    s->computeElementContributions_Acoustic_Flat(data);  // warmup
    FENCE
    for (auto _ : state)
      state.SetIterationTime(timeDevice([&] { s->computeElementContributions_Acoustic_Flat(data); }));
  }

  setThroughput(state, this->n_dof);
  this->setLabel(state);
}

// ==============================================================================
// BENCHMARK 2 : Teams (makutu)
// ==============================================================================
BENCHMARK_TEMPLATE_METHOD_F(AcousticKernelFixture, ElementContrib_Teams)(benchmark::State& state) {
  auto model = this->createModel();
  auto solver =
      solver_factory::createSolver(methodType::kSem, this->implem_, meshType::kUnstruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kAcoustic, this->order);
  solver->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);

  BenchmarkArrays arrays(this->n_rhs, this->n_time_steps, this->n_dof, model->getNumberOfPointsPerElement());
  auto wavefield = WavefieldAcoustic(arrays.pnGlobalPrev, arrays.pnGlobalCurr);
  auto rhs = RhsAcoustic(arrays.rhsTerm, arrays.rhsElement, arrays.rhsWeights);
  SEMsolverDataAcoustic data(wavefield, rhs);

  using Integral = typename IntegralTypeSelector<this->order, IntegralType::MAKUTU>::type;
  using Model = ModelUnstruct<float, int>;
  using SolverOnNodes = solver::fe::SEMsolver<this->order, Integral, Model, true, solver::fe::physicType::kAcoustic>;
  using SolverOnElems = solver::fe::SEMsolver<this->order, Integral, Model, false, solver::fe::physicType::kAcoustic>;

  if (this->isModelOnNodes_) {
    auto* s = dynamic_cast<SolverOnNodes*>(solver.get());
    s->computeElementContributions_Acoustic_Teams(data);  // warmup
    FENCE
    for (auto _ : state)
      state.SetIterationTime(timeDevice([&] { s->computeElementContributions_Acoustic_Teams(data); }));
  } else {
    auto* s = dynamic_cast<SolverOnElems*>(solver.get());
    s->computeElementContributions_Acoustic_Teams(data);  // warmup
    FENCE
    for (auto _ : state)
      state.SetIterationTime(timeDevice([&] { s->computeElementContributions_Acoustic_Teams(data); }));
  }

  setThroughput(state, this->n_dof);
  this->setLabel(state);
}

// ==============================================================================
// BENCHMARK 3 : Tensorial GEMM  (+ correctness gate vs Flat)
// ==============================================================================
BENCHMARK_TEMPLATE_METHOD_F(AcousticKernelFixture, ElementContrib_Tensorial)(benchmark::State& state) {
  auto model = this->createModel();

  using GemmIntegral = typename Qk_Hexahedron_Tensorial_Selector<this->order>::type;
  using RefIntegral = typename IntegralTypeSelector<this->order, IntegralType::MAKUTU>::type;
  using Model = ModelUnstruct<float, int>;
  using GemmOnNodes = solver::fe::SEMsolver<this->order, GemmIntegral, Model, true, solver::fe::physicType::kAcoustic>;
  using GemmOnElems = solver::fe::SEMsolver<this->order, GemmIntegral, Model, false, solver::fe::physicType::kAcoustic>;
  using RefOnNodes = solver::fe::SEMsolver<this->order, RefIntegral, Model, true, solver::fe::physicType::kAcoustic>;
  using RefOnElems = solver::fe::SEMsolver<this->order, RefIntegral, Model, false, solver::fe::physicType::kAcoustic>;

  BenchmarkArrays arrays(this->n_rhs, this->n_time_steps, this->n_dof, model->getNumberOfPointsPerElement());
  auto wavefield = WavefieldAcoustic(arrays.pnGlobalPrev, arrays.pnGlobalCurr);
  auto rhs = RhsAcoustic(arrays.rhsTerm, arrays.rhsElement, arrays.rhsWeights);
  SEMsolverDataAcoustic data(wavefield, rhs);

  // Reference makutu Flat solver, used once for the correctness gate.
  auto ref =
      solver_factory::createSolver(methodType::kSem, implemType::kMakutu, meshType::kUnstruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kAcoustic, this->order);
  ref->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);

  constexpr double kRelTol = 1e-4;

  if (this->isModelOnNodes_) {
    GemmOnNodes my_solver;
    my_solver.computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);
    auto* r = dynamic_cast<RefOnNodes*>(ref.get());

    // --- correctness gate (once, untimed): GEMM must match Flat ---
    r->computeElementContributions_Acoustic_Flat(data);
    my_solver.computeElementContributions_Acoustic_Gemm(data);
    Kokkos::fence();
    const double rel = maxRelDiffForce(*r, my_solver, 1);
    state.counters["relDiffVsFlat"] = rel;
    if (!(rel < kRelTol)) {
      state.SkipWithError("GEMM result diverges from Flat (see relDiffVsFlat)");
      return;
    }

    my_solver.computeElementContributions_Acoustic_Gemm(data);  // warmup
    FENCE
    for (auto _ : state)
      state.SetIterationTime(timeDevice([&] { my_solver.computeElementContributions_Acoustic_Gemm(data); }));
  } else {
    GemmOnElems my_solver;
    my_solver.computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);
    auto* r = dynamic_cast<RefOnElems*>(ref.get());

    r->computeElementContributions_Acoustic_Flat(data);
    my_solver.computeElementContributions_Acoustic_Gemm(data);
    Kokkos::fence();
    const double rel = maxRelDiffForce(*r, my_solver, 1);
    state.counters["relDiffVsFlat"] = rel;
    if (!(rel < kRelTol)) {
      state.SkipWithError("GEMM result diverges from Flat (see relDiffVsFlat)");
      return;
    }

    my_solver.computeElementContributions_Acoustic_Gemm(data);  // warmup
    FENCE
    for (auto _ : state)
      state.SetIterationTime(timeDevice([&] { my_solver.computeElementContributions_Acoustic_Gemm(data); }));
  }

  setThroughput(state, this->n_dof);
  this->setLabel(state);
}

BENCHMARK_FOR_ALL_ORDERS(AcousticKernelFixture, ElementContrib_Flat,
                         AcousticKernelConfig,
                             ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})
                             ->UseManualTime()
                             ->Unit(benchmark::kMicrosecond))

BENCHMARK_FOR_ALL_ORDERS(AcousticKernelFixture, ElementContrib_Teams,
                         AcousticKernelConfig,
                             ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})
                             ->UseManualTime()
                             ->Unit(benchmark::kMicrosecond))

BENCHMARK_FOR_ALL_ORDERS(AcousticKernelFixture, ElementContrib_Tensorial,
                         AcousticKernelConfig,
                             ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})
                             ->UseManualTime()
                             ->Unit(benchmark::kMicrosecond))

}  // namespace bench
}  // namespace model

int main(int argc, char** argv) { return runBenchmarks(argc, argv); }
