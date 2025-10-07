#pragma once

/**
 * @brief Runs Google Benchmark benchmarks with optional Kokkos initialization.
 *
 * This function serves as the main entry point for running benchmarks. It handles
 * initialization and finalization of both Google Benchmark and Kokkos (if enabled).
 * When USE_KOKKOS is defined, Kokkos is initialized before any benchmarks run and
 * finalized after all benchmarks complete.
 *
 * @param argc The number of command-line arguments passed to the program.
 * @param argv The array of command-line argument strings.
 *
 * @return 0 if benchmarks ran successfully, 1 if unrecognized arguments were provided.
 *
 * @note If USE_KOKKOS is defined, Kokkos::initialize() and Kokkos::finalize() are
 *       called exactly once before and after all benchmarks, respectively.
 * @note This function ensures proper cleanup (Kokkos finalization) even when
 *       unrecognized arguments are detected.
 */
static int runBenchmarks(int argc, char** argv)
{
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