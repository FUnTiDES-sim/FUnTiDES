#include "Integrals.h"
#include "bench_main.h"
#include "bench_solver_fe_cartesian_unstruct_common.h"
#include "model_unstruct.h"
#include "rhs_acoustic.h"
#include "sem_enums.h"
#include "sem_solver.h"
#include "wavefield_acoustic.h"

namespace model {
namespace bench {

template <int Order>
using AcousticKernelConfig = BuilderConfig<Order, false>;

// Alias pour éviter les conflits avec ton autre fichier de bench
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

// ==============================================================================
// BENCHMARK 1 : Flat
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
    auto* solver_c = dynamic_cast<SolverOnNodes*>(solver.get());
    for (auto _ : state) {
      solver_c->computeElementContributions_Acoustic_Flat(data);
      FENCE
    }
  } else {
    auto* solver_c = dynamic_cast<SolverOnElems*>(solver.get());
    for (auto _ : state) {
      solver_c->computeElementContributions_Acoustic_Flat(data);
      FENCE
    }
  }

  this->setLabel(state);
}

// ==============================================================================
// BENCHMARK 2 : Teams
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
    auto* solver_c = dynamic_cast<SolverOnNodes*>(solver.get());
    for (auto _ : state) {
      solver_c->computeElementContributions_Acoustic_Teams(data);
      FENCE
    }
  } else {
    auto* solver_c = dynamic_cast<SolverOnElems*>(solver.get());
    for (auto _ : state) {
      solver_c->computeElementContributions_Acoustic_Teams(data);
      FENCE
    }
  }

  this->setLabel(state);
}

BENCHMARK_FOR_ALL_ORDERS(AcousticKernelFixture, ElementContrib_Flat,
                         AcousticKernelConfig,
                             ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})
                             ->Unit(benchmark::kMicrosecond))

BENCHMARK_FOR_ALL_ORDERS(AcousticKernelFixture, ElementContrib_Teams,
                         AcousticKernelConfig,
                             ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})
                             ->Unit(benchmark::kMicrosecond))

}  // namespace bench
}  // namespace model

int main(int argc, char** argv) { return runBenchmarks(argc, argv); }
