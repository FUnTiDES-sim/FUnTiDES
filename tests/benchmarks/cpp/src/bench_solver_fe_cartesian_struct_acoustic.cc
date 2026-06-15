#include "bench_main.h"
#include "bench_solver_fe_cartesian_struct_common.h"
#include "rhs_acoustic.h"
#include "wavefield_acoustic.h"

namespace model {
namespace bench {

template <int Order>
using AcousticConfig = BuilderConfig<Order, false>;

// Alias to avoid any benchmark naming conflict
template <typename T>
class AcousticFixture : public SolverStructFixture<T> {};

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

BENCHMARK_TEMPLATE_METHOD_F(AcousticFixture, FEInit)
(benchmark::State& state) {
  auto model = this->createModel();

  auto solver =
      solver_factory::createSolver(methodType::kSem, this->implem_, meshType::kStruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kAcoustic, this->order);

  for (auto _ : state) {
    solver->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);
  }

  this->setLabel(state);
}

BENCHMARK_TEMPLATE_METHOD_F(AcousticFixture, OneStep)
(benchmark::State& state) {
  auto model = this->createModel();

  auto solver =
      solver_factory::createSolver(methodType::kSem, this->implem_, meshType::kStruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kAcoustic, this->order);

  solver->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);

  BenchmarkArrays arrays(this->n_rhs, this->n_time_steps, this->n_dof, model->getNumberOfPointsPerElement());

  arrays.rhsElement(0) = this->ex / 2 + this->ey / 2 * this->ex + this->ez / 2 * this->ey * this->ex;
  arrays.rhsElement(1) = this->ex / 3 + this->ey / 2 * this->ex + this->ez / 2 * this->ey * this->ex;

  SolverUtils myUtils;
  float const tpeak = 1.0f / this->f0;
  std::vector<float> sourceTerm = myUtils.computeSourceTerm(this->n_time_steps, this->dt, this->f0, 2, tpeak);
  for (int j = 0; j < this->n_time_steps; j++) {
    arrays.rhsTerm(0, j) = sourceTerm[j];
  }

  auto wavefield = WavefieldAcoustic(arrays.pnGlobalPrev, arrays.pnGlobalCurr);
  auto rhs = RhsAcoustic(arrays.rhsTerm, arrays.rhsElement, arrays.rhsWeights);
  SEMsolverDataAcoustic data(wavefield, rhs);

  for (auto _ : state) {
    solver->computeForces(this->dt, this->time_sample, data);
    solver->updateSolutionForward(this->dt, data);
  }

  this->setLabel(state);
}

BENCHMARK_FOR_ALL_ORDERS(
    AcousticFixture, FEInit,
    AcousticConfig, ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})->Unit(benchmark::kMillisecond))
BENCHMARK_FOR_ALL_ORDERS(
    AcousticFixture, OneStep,
    AcousticConfig, ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})->Unit(benchmark::kMillisecond))

}  // namespace bench
}  // namespace model

int main(int argc, char** argv) { return runBenchmarks(argc, argv); }
