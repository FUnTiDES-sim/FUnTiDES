#ifndef TESTS_BENCHMARKS_CPP_INCLUDE_BENCH_SOLVER_FE_CARTESIAN_STRUCT_COMMON_H_
#define TESTS_BENCHMARKS_CPP_INCLUDE_BENCH_SOLVER_FE_CARTESIAN_STRUCT_COMMON_H_

#include <benchmark/benchmark.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "bench_macros.h"
#include "cartesian_struct_builder.h"
#include "data_type.h"
#include "model.h"
#include "sem_solver.h"
#include "solver_factory.h"
#include "utils.h"

using namespace solver::fe;
using namespace utils::enums;

namespace model {
namespace bench {

// Template config for the fixture (factorized)
template <int Order, bool IsElastic>
struct BuilderConfig {
  using Builder = CartesianStructBuilder<float, int, Order>;
  static constexpr int order = Order;
  static constexpr bool is_elastic = IsElastic;
};

// Base Template fixture for the benchmarks
template <typename T>
class SolverStructFixture : public benchmark::Fixture {
 protected:
  // Avoid hiding warnings for overloaded virtual functions
  using benchmark::Fixture::SetUp;
  using benchmark::Fixture::TearDown;

  // domain decomposition
  static constexpr int rank = 0;
  static constexpr int size = 1;
  static constexpr float origin = 0.0f;
  float local_l = 2000.0f;

  // model
  static constexpr int ex = 100;
  static constexpr int ey = 100;
  static constexpr int ez = 100;
  static constexpr float domain_size = 2000.0f;
  static constexpr int order = T::order;
  static constexpr bool is_elastic = T::is_elastic;
  static constexpr int n_dof = (ex * order + 1) * (ey * order + 1) * (ez * order + 1);
  bool isModelOnNodes_;

  // sponge
  inline static constexpr std::array<float, 3> sponge_size = {200.0f, 200.0f, 200.0f};
  inline static constexpr bool surface_sponge = false;
  inline static constexpr float taper_delta = 100.0f;

  // solver
  static constexpr int n_rhs = 2;
  static constexpr float dt = 0.001f;
  static constexpr int time_sample = 1;
  static constexpr int n_time_steps = 1500;
  static constexpr float f0 = 5.0f;
  implemType implem_;

  void SetUp(const ::benchmark::State& state) override {
    isModelOnNodes_ = state.range(0);
    implem_ = static_cast<implemType>(state.range(1));
    local_l = domain_size;
  }

  std::shared_ptr<model::ModelApi<float, int>> createModel() {
    float hx = domain_size / ex;
    float hy = domain_size / ey;
    float hz = domain_size / ez;

    // Note: Builder defaults origins to 0.0, which matches our serial setup
    typename T::Builder builder(ex, hx, ey, hy, ez, hz, isModelOnNodes_, is_elastic);
    return builder.getModel(true);
  }

  void setLabel(benchmark::State& state) const {
    state.SetLabel("Order=" + std::to_string(order) + " OnNodes=" + to_string(isModelOnNodes_) +
                   " Implem=" + to_string(implem_) + " IsElastic=" + std::to_string(is_elastic));
  }
};

}  // namespace bench
}  // namespace model

#endif  // TESTS_BENCHMARKS_CPP_INCLUDE_BENCH_SOLVER_FE_CARTESIAN_STRUCT_COMMON_H_
