#pragma once

#include <benchmark/benchmark.h>

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

/**
 * @brief Parses the --kokkos-threads argument from command line arguments.
 * @param argc Argument count
 * @param argv Argument vector
 * @return Number of Kokkos threads specified
 * @throws std::invalid_argument if argument is missing, malformed, or negative
 */
static int parseKokkosThreads(int argc, char** argv)
{
  const char prefix[] = "--kokkos-threads=";
  const std::size_t len = sizeof(prefix) - 1;
  for (int i = 1; i < argc; ++i)
  {
    const char* arg = argv[i];
    if (std::strncmp(arg, prefix, len) == 0)
    {
      const char* value = arg + len;
      if (*value == '\0')
        throw std::invalid_argument("--kokkos-threads requires a value");
      char* end = nullptr;
      long v = std::strtol(value, &end, 10);
      if (*end != '\0' || v < 0)
        throw std::invalid_argument(
            "--kokkos-threads must be a non-negative integer");
      return static_cast<int>(v);
    }
  }
  throw std::invalid_argument(
      "--kokkos-threads must be set to > 0, e.g. --kokkos-threads=4");
}

/**
 * @brief Runs Google Benchmark benchmarks with optional Kokkos initialization.
 *
 * This function serves as the main entry point for running benchmarks. It
 * handles initialization and finalization of both Google Benchmark and Kokkos
 * (if enabled). When USE_KOKKOS is defined, Kokkos is initialized before any
 * benchmarks run and finalized after all benchmarks complete.
 *
 * @param argc The number of command-line arguments passed to the program.
 * @param argv The array of command-line argument strings. Must include
 * --kokkos-threads if USE_KOKKOS is defined.
 *
 * @return 0 if benchmarks ran successfully, 1 if unrecognized arguments were
 * provided.
 *
 * @note If USE_KOKKOS is defined, Kokkos::initialize() and Kokkos::finalize()
 * are called exactly once before and after all benchmarks, respectively.
 * @note This function ensures proper cleanup (Kokkos finalization) even when
 *       unrecognized arguments are detected.
 */
static int runBenchmarks(int argc, char** argv)
{
#ifdef USE_KOKKOS
  int nThreads = parseKokkosThreads(argc, argv);
  Kokkos::InitializationSettings kkSettings;
  kkSettings.set_num_threads(nThreads);
  Kokkos::initialize(kkSettings);
#endif

  benchmark::Initialize(&argc, argv);
#ifdef USE_KOKKOS
  benchmark::AddCustomContext("kokkos_threads", std::to_string(nThreads));
#endif
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();

#ifdef USE_KOKKOS
  // Finalize Kokkos ONCE after all benchmarks
  Kokkos::finalize();
#endif

  return 0;
}