#pragma once

/**
 * @brief Macro to instantiate benchmark templates for multiple template order
 * values.
 *
 * This macro generates four benchmark instantiations for orders 1 through 3,
 * where each instantiation uses the specified fixture, name, and type template
 * with a configurable configuration suffix.
 *
 * @param FIXTURE The benchmark fixture class to use for the instantiations
 * @param NAME The name of the benchmark to instantiate
 * @param TYPE The template type that accepts an integer template parameter
 * (order)
 * @param CONFIG Additional configuration or method calls to append to each
 * instantiation
 *
 * @note The CONFIG parameter should include any leading dots or arrows for
 * method chaining (e.g., "->Unit(benchmark::kMillisecond)")
 *
 * @note Known as the lizard macro ;)
 *
 * Example usage:
 * @code
 * BENCHMARK_FOR_ALL_ORDERS(MyFixture, MyBenchmark, MyType,
 * ->Unit(benchmark::kMicrosecond))
 * @endcode
 */
#define BENCHMARK_FOR_ALL_ORDERS(FIXTURE, NAME, TYPE, CONFIG)      \
  BENCHMARK_TEMPLATE_INSTANTIATE_F(FIXTURE, NAME, TYPE<1>) CONFIG; \
  BENCHMARK_TEMPLATE_INSTANTIATE_F(FIXTURE, NAME, TYPE<2>) CONFIG; \
  BENCHMARK_TEMPLATE_INSTANTIATE_F(FIXTURE, NAME, TYPE<3>) CONFIG;
