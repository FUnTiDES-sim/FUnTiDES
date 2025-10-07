#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <memory>

#include "utils.h"
#include "cartesian_struct_builder.h"
#include "data_type.h"
#include "model.h"
#include "sem_solver.h"
#include "solver_base.h"
#include "solver_factory.h"

// Constants
constexpr float PI = 3.14157f;

// Structure to hold test data configuration
struct StructData
{
  int ex, ey, ez;
  float domain_size;
  float hx, hy, hz;
  int order;
  int nx, ny, nz;
  int n_dof;

  StructData(int order_)
      : ex(100),
        ey(100),
        ez(100),
        domain_size(2000.0f),
        order(order_)
  {
    hx = domain_size / ex;
    hy = domain_size / ey;
    hz = domain_size / ez;
    nx = ex * order + 1;
    ny = ey * order + 1;
    nz = ez * order + 1;
    n_dof = nx * ny * nz;
  }
};

// Benchmark fixture
class SolverStructFixture : public benchmark::Fixture
{
 protected:
  static constexpr int ex = 100;
  static constexpr int ey = 100;
  static constexpr int ez = 100;
  static constexpr float domain_size = 2000.0f;

  std::shared_ptr<model::ModelApi<float, int>> createModel(
      int order, bool isModelOnNodes)
  {
    float hx = domain_size / ex;
    float hy = domain_size / ey;
    float hz = domain_size / ez;

    switch (order)
    {
      case 1:
      {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, hx, ey, hy, ez, hz, isModelOnNodes);
        return builder.getModel();
      }
      case 2:
      {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, hx, ey, hy, ez, hz, isModelOnNodes);
        return builder.getModel();
      }
      case 3:
      {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, hx, ey, hy, ez, hz, isModelOnNodes);
        return builder.getModel();
      }
      case 4:
      {
        model::CartesianStructBuilder<float, int, 4> builder(
            ex, hx, ey, hy, ez, hz, isModelOnNodes);
        return builder.getModel();
      }
      default:
        return nullptr;
    }
  }
};

// Structure to hold allocated arrays for benchmarks
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
    rhsWeights = allocateArray2D<arrayReal>(n_rhs, nb_points_per_element, "rhsWeights");
    pnGlobal = allocateArray2D<arrayReal>(n_dof, 2, "pnGlobal");
    rhsLocation = allocateArray2D<arrayReal>(1, 3, "rhsLocation");
    
    FENCE
  }
};

// Benchmark for compute_fe_init using fixture
BENCHMARK_DEFINE_F(SolverStructFixture, FEInit)(benchmark::State &state)
{
  int order = state.range(0);
  bool isModelOnNodes = state.range(1);
  SolverFactory::implemType implem = SolverFactory::implemType::MAKUTU;

  // Create model
  auto model = createModel(order, isModelOnNodes);

  // Sponge parameters
  std::array<float, 3> sponge_size = {200.0f, 200.0f, 200.0f};
  bool surface_sponge = false;
  float taper_delta = 100.0f;

  // Create solver
  auto solver = SolverFactory::createSolver(
      SolverFactory::methodType::SEM, implem,
      SolverFactory::meshType::Struct, order);

  for (auto _ : state)
  {
    // Benchmark the initialization
    solver->computeFEInit(*model, sponge_size, surface_sponge, taper_delta);
  }

  // Report additional metrics
  state.SetLabel("Order=" + std::to_string(order) +
                 " OnNodes=" + std::to_string(isModelOnNodes));
}

// Benchmark for compute_one_step using fixture
BENCHMARK_DEFINE_F(SolverStructFixture, OneStep)(benchmark::State &state)
{
  int order = state.range(0);
  bool isModelOnNodes = state.range(1);
  SolverFactory::implemType implem = SolverFactory::implemType::MAKUTU;

  StructData sd(order);

  // Parameters
  int n_rhs = 2;
  float dt = 0.001f;
  int time_sample = 1;
  int n_time_steps = 1500;
  float f0 = 5.0f;

  // Create model
  auto model = createModel(order, isModelOnNodes);

  // Sponge parameters
  std::array<float, 3> sponge_size = {200.0f, 200.0f, 200.0f};
  bool surface_sponge = false;
  float taper_delta = 100.0f;

  // Create and initialize solver
  auto solver = SolverFactory::createSolver(
      SolverFactory::methodType::SEM, implem,
      SolverFactory::meshType::Struct, order);

  solver->computeFEInit(*model, sponge_size, surface_sponge, taper_delta);

  BenchmarkArrays arrays(n_rhs, n_time_steps, sd.n_dof, 
                         model->getNumberOfPointsPerElement());

  arrays.rhsElement(0) = ex / 2 + ey / 2 * ex + ez / 2 * ey * ex;
  arrays.rhsElement(1) = ex / 3 + ey / 2 * ex + ez / 2 * ey * ex;

  SolverUtils myUtils;
  std::vector<float> sourceTerm = myUtils.computeSourceTerm(n_time_steps, dt, f0, 2);
  for (int j = 0; j < n_time_steps; j++)
  {
    arrays.rhsTerm(0, j) = sourceTerm[j];
  }

  // Create solver data
  SEMsolverData data(0, 1, arrays.rhsTerm, arrays.pnGlobal, arrays.rhsElement, arrays.rhsWeights);

  for (auto _ : state)
  {
    // Benchmark one time step
    solver->computeOneStep(dt, time_sample, data);
  }

}

// Register benchmarks with different parameters
// FEInit benchmarks
//BENCHMARK_REGISTER_F(SolverStructFixture, FEInit)->Args({1, 1})->ThreadRange(1, 8)->Unit(benchmark::kMillisecond);
//BENCHMARK_REGISTER_F(SolverStructFixture, FEInit)->Args({1, 0})->Unit(benchmark::kMillisecond);
//BENCHMARK_REGISTER_F(SolverStructFixture, FEInit)->Args({2, 1})->Unit(benchmark::kMillisecond);
//BENCHMARK_REGISTER_F(SolverStructFixture, FEInit)->Args({2, 0})->Unit(benchmark::kMillisecond);
//BENCHMARK_REGISTER_F(SolverStructFixture, FEInit)->Args({3, 1})->Unit(benchmark::kMillisecond);
//BENCHMARK_REGISTER_F(SolverStructFixture, FEInit)->Args({3, 0})->Unit(benchmark::kMillisecond);

// OneStep benchmarks
BENCHMARK_REGISTER_F(SolverStructFixture, OneStep)->Args({1, 1})->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(SolverStructFixture, OneStep)->Args({1, 0})->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(SolverStructFixture, OneStep)->Args({2, 1})->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(SolverStructFixture, OneStep)->Args({2, 0})->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(SolverStructFixture, OneStep)->Args({3, 1})->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(SolverStructFixture, OneStep)->Args({3, 0})->Unit(benchmark::kMillisecond);

int main(int argc, char** argv) {
#ifdef USE_KOKKOS
  // Initialize Kokkos ONCE before any benchmarks
  Kokkos::initialize(argc, argv);
#endif
  
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
#ifdef USE_KOKKOS
    Kokkos::finalize();
#endif
    return 1;
  }
  
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  
#ifdef USE_KOKKOS
  // Finalize Kokkos ONCE after all benchmarks
  Kokkos::finalize();
#endif
  
  return 0;
}