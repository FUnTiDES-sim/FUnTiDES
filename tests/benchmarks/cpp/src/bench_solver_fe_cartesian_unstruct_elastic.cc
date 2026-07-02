#include "bench_main.h"
#include "bench_solver_fe_cartesian_unstruct_common.h"
#include "rhs_elastic.h"
#include "wavefield_elastic.h"

namespace model {
namespace bench {

template <int Order>
using ElasticConfig = BuilderConfig<Order, true>;

// Alias to avoid any benchmark naming conflict
template <typename T>
class ElasticFixture : public SolverUnstructFixture<T> {};

struct BenchmarkArrays {
  arrayReal rhsTermx;
  arrayReal rhsTermy;
  arrayReal rhsTermz;
  vectorInt rhsElement;
  arrayReal rhsWeights;
  vectorReal uxnGlobalPrev;
  vectorReal uynGlobalPrev;
  vectorReal uznGlobalPrev;
  vectorReal uxnGlobalCurr;
  vectorReal uynGlobalCurr;
  vectorReal uznGlobalCurr;
  arrayReal rhsLocation;

  BenchmarkArrays(int n_rhs, int n_time_steps, int n_dof, int nb_points_per_element) {
    rhsTermx = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTermx");
    rhsTermy = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTermy");
    rhsTermz = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTermz");
    rhsElement = allocateVector<vectorInt>(n_rhs, "rhsElement");
    rhsWeights = allocateArray2D<arrayReal>(n_rhs, nb_points_per_element, "rhsWeights");
    uxnGlobalPrev = allocateVector<vectorReal>(n_dof, "uxnGlobalPrev");
    uynGlobalPrev = allocateVector<vectorReal>(n_dof, "uynGlobalPrev");
    uznGlobalPrev = allocateVector<vectorReal>(n_dof, "uznGlobalPrev");
    uxnGlobalCurr = allocateVector<vectorReal>(n_dof, "uxnGlobalCurr");
    uynGlobalCurr = allocateVector<vectorReal>(n_dof, "uynGlobalCurr");
    uznGlobalCurr = allocateVector<vectorReal>(n_dof, "uznGlobalCurr");
    rhsLocation = allocateArray2D<arrayReal>(1, 3, "rhsLocation");

    FENCE
  }
};

BENCHMARK_TEMPLATE_METHOD_F(ElasticFixture, FEInit)
(benchmark::State& state) {
  auto model = this->createModel();
  auto solver =
      solver_factory::createSolver(methodType::kSem, this->implem_, meshType::kUnstruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kElastic, this->order);

  for (auto _ : state) {
    solver->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);
  }

  this->setLabel(state);
}

BENCHMARK_TEMPLATE_METHOD_F(ElasticFixture, OneStep)
(benchmark::State& state) {
  auto model = this->createModel();
  auto solver =
      solver_factory::createSolver(methodType::kSem, this->implem_, meshType::kUnstruct,
                                   this->isModelOnNodes_ ? modelLocationType::kOnNodes : modelLocationType::kOnElements,
                                   physicType::kElastic, this->order);

  solver->computeFEInit(*model, this->sponge_size, this->surface_sponge, this->taper_delta);

  BenchmarkArrays arrays(this->n_rhs, this->n_time_steps, this->n_dof, model->getNumberOfPointsPerElement());

  arrays.rhsElement(0) = this->ex / 2 + this->ey / 2 * this->ex + this->ez / 2 * this->ey * this->ex;
  arrays.rhsElement(1) = this->ex / 3 + this->ey / 2 * this->ex + this->ez / 2 * this->ey * this->ex;

  SolverUtils myUtils;
  float const tpeak = 1.0f / this->f0;
  std::vector<float> sourceTerm = myUtils.computeSourceTerm(this->n_time_steps, this->dt, this->f0, 2, tpeak);
  for (int j = 0; j < this->n_time_steps; j++) {
    arrays.rhsTermx(0, j) = sourceTerm[j];
    arrays.rhsTermy(0, j) = sourceTerm[j];
    arrays.rhsTermz(0, j) = sourceTerm[j];
  }

  auto wavefield = WavefieldElastic(arrays.uxnGlobalPrev, arrays.uynGlobalPrev, arrays.uznGlobalPrev,
                                    arrays.uxnGlobalCurr, arrays.uynGlobalCurr, arrays.uznGlobalCurr);
  auto rhs = RhsElastic(arrays.rhsTermx, arrays.rhsTermy, arrays.rhsTermz, arrays.rhsElement, arrays.rhsWeights);
  SEMsolverDataElastic data(wavefield, rhs);

  for (auto _ : state) {
    solver->computeForces(this->dt, this->time_sample, data);
    solver->updateSolutionForward(this->dt, data);
  }

  this->setLabel(state);
}

BENCHMARK_FOR_ALL_ORDERS(
    ElasticFixture, FEInit,
    ElasticConfig, ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})->Unit(benchmark::kMillisecond))
BENCHMARK_FOR_ALL_ORDERS(
    ElasticFixture, OneStep,
    ElasticConfig, ->ArgsProduct({{0, 1}, {static_cast<int64_t>(implemType::kMakutu)}})->Unit(benchmark::kMillisecond))

}  // namespace bench
}  // namespace model

int main(int argc, char** argv) { return runBenchmarks(argc, argv); }
