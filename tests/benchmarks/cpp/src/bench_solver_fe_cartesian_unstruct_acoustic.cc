#include <benchmark/benchmark.h>

#include <array>
#include <memory>

#include "bench_macros.h"
#include "bench_main.h"
#include "cartesian_unstruct_builder.h"
#include "data_type.h"
#include "model.h"
#include "sem_solver_acoustic.h"
#include "solver_factory.h"
#include "utils.h"

namespace model
{
namespace bench
{

// Template config for the fixture
template <int Order>
struct BuilderConfig
{
  using Builder = CartesianUnstructBuilder<float, int>;
  using BuilderParams = CartesianParams<float, int>;
  static constexpr int order = Order;
};

// Template fixture for the benchmarks
template <typename T>
class SolverUnstructFixture : public benchmark::Fixture
{
 protected:
  // model
  static constexpr int ex = 100;
  static constexpr int ey = 100;
  static constexpr int ez = 100;
  static constexpr float lx = 2000.0f;
  static constexpr float ly = 2000.0f;
  static constexpr float lz = 2000.0f;
  static constexpr int order = T::order;
  static constexpr int n_dof =
      (ex * order + 1) * (ey * order + 1) * (ez * order + 1);
  bool isModelOnNodes_;

  // sponge
  inline static constexpr std::array<float, 3> sponge_size = {200.0f, 200.0f,
                                                              200.0f};
  inline static constexpr bool surface_sponge = false;
  inline static constexpr float taper_delta = 100.0f;

  // solver
  static constexpr int n_rhs = 2;
  static constexpr float dt = 0.001f;
  static constexpr int time_sample = 1;
  static constexpr int n_time_steps = 1500;
  static constexpr float f0 = 5.0f;
  SolverFactory::implemType implem_;

  void SetUp(const ::benchmark::State& state) override
  {
    isModelOnNodes_ = state.range(0);
    implem_ = static_cast<SolverFactory::implemType>(state.range(1));
  }

  std::shared_ptr<model::ModelApi<float, int>> createModel()
  {
    typename T::BuilderParams params(order, ex, ey, ez, lx, ly, lz,
                                     isModelOnNodes_, false);
    typename T::Builder builder(params);
    return builder.getModel();
  }

  void setLabel(benchmark::State& state) const
  {
    state.SetLabel("Order=" + std::to_string(order) +
                   " OnNodes=" + std::to_string(isModelOnNodes_) +
                   " Implem=" + std::to_string(implem_) +
                   " IsElastic=" + std::to_string(false));
  }
};

// Unstructure to hold allocated arrays for benchmarks
struct BenchmarkArrays
{
  arrayReal rhsTerm;
  vectorInt rhsElement;
  arrayReal rhsWeights;
  arrayReal pnGlobal;
  arrayReal rhsLocation;

  BenchmarkArrays(int n_rhs, int n_time_steps, int n_dof,
                  int nb_points_per_element)
  {
    rhsTerm = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTerm");
    rhsElement = allocateVector<vectorInt>(n_rhs, "rhsElement");
    rhsWeights =
        allocateArray2D<arrayReal>(n_rhs, nb_points_per_element, "rhsWeights");
    pnGlobal = allocateArray2D<arrayReal>(n_dof, 2, "pnGlobal");
    rhsLocation = allocateArray2D<arrayReal>(1, 3, "rhsLocation");

    FENCE
  }
};

BENCHMARK_TEMPLATE_METHOD_F(SolverUnstructFixture, FEInit)
(benchmark::State& state)
{
  // Prepare
  auto model = this->createModel();

  auto solver = SolverFactory::createSolver(
      SolverFactory::methodType::SEM, this->implem_,
      SolverFactory::meshType::Unstruct,
      this->isModelOnNodes_ ? SolverFactory::modelLocationType::OnNodes
                            : SolverFactory::modelLocationType::OnElements,
      SolverFactory::physicType::Acoustic, this->order);

  // Bench
  for (auto _ : state)
  {
    solver->computeFEInit(*model, this->sponge_size, this->surface_sponge,
                          this->taper_delta);
  }

  // Label
  this->setLabel(state);
}

BENCHMARK_TEMPLATE_METHOD_F(SolverUnstructFixture, OneStep)
(benchmark::State& state)
{
  // Prepare
  auto model = this->createModel();

  auto solver = SolverFactory::createSolver(
      SolverFactory::methodType::SEM, this->implem_,
      SolverFactory::meshType::Unstruct,
      this->isModelOnNodes_ ? SolverFactory::modelLocationType::OnNodes
                            : SolverFactory::modelLocationType::OnElements,
      SolverFactory::physicType::Acoustic, this->order);

  solver->computeFEInit(*model, this->sponge_size, this->surface_sponge,
                        this->taper_delta);

  BenchmarkArrays arrays(this->n_rhs, this->n_time_steps, this->n_dof,
                         model->getNumberOfPointsPerElement());
  // sources at the center of the domain
  arrays.rhsElement(0) = this->ex / 2 + this->ey / 2 * this->ex +
                         this->ez / 2 * this->ey * this->ex;
  arrays.rhsElement(1) = this->ex / 3 + this->ey / 2 * this->ex +
                         this->ez / 2 * this->ey * this->ex;

  // ricker wavelet
  SolverUtils myUtils;
  std::vector<float> sourceTerm =
      myUtils.computeSourceTerm(this->n_time_steps, this->dt, this->f0, 2);
  for (int j = 0; j < this->n_time_steps; j++)
  {
    arrays.rhsTerm(0, j) = sourceTerm[j];
  }

  SEMsolverDataAcoustic data(0, 1, arrays.rhsTerm, arrays.pnGlobal,
                             arrays.rhsElement, arrays.rhsWeights);

  // Bench
  for (auto _ : state)
  {
    solver->computeOneStep(this->dt, this->time_sample, data);
  }

  // Label
  this->setLabel(state);
}

// Instantiate for all order/isModelOnNodes/implemType combinations
// TODO add SolverFactory::implemType::SHIVA when reactivated in compilation
BENCHMARK_FOR_ALL_ORDERS(
    SolverUnstructFixture, FEInit,
    BuilderConfig,
        ->ArgsProduct({{0, 1}, {SolverFactory::implemType::MAKUTU}})
        ->Unit(benchmark::kMillisecond))
BENCHMARK_FOR_ALL_ORDERS(
    SolverUnstructFixture, OneStep,
    BuilderConfig,
        ->ArgsProduct({{0, 1}, {SolverFactory::implemType::MAKUTU}})
        ->Unit(benchmark::kMillisecond))

}  // namespace bench
}  // namespace model

int main(int argc, char** argv) { return runBenchmarks(argc, argv); }