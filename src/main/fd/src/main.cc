//************************************************************************
// Finite Difference Time Domain (FDTD) Acoustic Simulation
// Version 0.0.1
//
// src/main/fd/main.cpp
//
// Main driver program for FDTD simulation
//
// This file provides the entry point for the FDTD acoustic wave
// propagation simulator. It handles command-line argument parsing,
// initializes the simulation environment (including optional Kokkos
// support), and orchestrates the complete simulation workflow.
//************************************************************************

#include <chrono>
#include <cstdlib>
#include <cxxopts.hpp>
#include <exception>
#include <iostream>
#include <memory>

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include "fd_options.h"
#include "fd_proxy.h"

namespace
{

using std::chrono::duration_cast;
using std::chrono::nanoseconds;
using std::chrono::system_clock;
using std::chrono::time_point;

// Global start time for total execution timing
time_point<system_clock> g_start_init_time;

#ifdef USE_KOKKOS
/**
 * @brief RAII wrapper for Kokkos initialization/finalization.
 */
class KokkosScope
{
 public:
  explicit KokkosScope(int argc, char* argv[])
  {
    Kokkos::initialize(argc, argv);
  }
  ~KokkosScope() { Kokkos::finalize(); }

  // Prevent copying and moving
  KokkosScope(const KokkosScope&) = delete;
  KokkosScope& operator=(const KokkosScope&) = delete;
  KokkosScope(KokkosScope&&) = delete;
  KokkosScope& operator=(KokkosScope&&) = delete;
};
#endif

/**
 * @brief Converts nanoseconds duration to seconds as a double.
 * @param duration Duration in nanoseconds
 * @return Duration in seconds
 */
inline double NanosecondsToSeconds(const nanoseconds& duration)
{
  return duration.count() / 1E9;
}

/**
 * @brief Prints a formatted timing report.
 * @param label Description of the timing measurement
 * @param duration Duration to report
 */
void PrintTiming(const std::string& label, const nanoseconds& duration)
{
  std::cout << label << ": " << NanosecondsToSeconds(duration) << " seconds."
            << std::endl;
}

/**
 * @brief Executes the complete FDTD simulation workflow.
 *
 * Initializes the FDTD simulation environment, runs the time-stepping
 * loop, and reports performance metrics.
 *
 * @param fd_sim Reference to the FdtdProxy simulation object
 */
void Compute(fdtd::FdtdProxy& fd_sim)
{
  // Initialize FDTD simulation
  std::cout << "Initializing FDTD simulation..." << std::endl;
  fd_sim.InitFdtd();
  std::cout << "FDTD initialization complete." << std::endl;

  // Start computation timer
  const auto start_run_time = system_clock::now();

  // Run simulation
  std::cout << "\nStarting FDTD computation..." << std::endl;
  fd_sim.Run();
  std::cout << "FDTD computation complete.\n" << std::endl;

  // Report timing information
  const auto init_duration =
      duration_cast<nanoseconds>(start_run_time - g_start_init_time);
  const auto compute_duration =
      duration_cast<nanoseconds>(system_clock::now() - start_run_time);

  std::cout << "=== Timing Summary ===" << std::endl;
  PrintTiming("Initialization time", init_duration);
  PrintTiming("Computation time   ", compute_duration);
  std::cout << "======================" << std::endl;
}

/**
 * @brief Configures environment variables for optimal Kokkos/OpenMP
 * performance.
 */
void ConfigureParallelEnvironment()
{
#ifdef USE_KOKKOS
  // Configure OpenMP thread binding for optimal performance
  setenv("OMP_PROC_BIND", "spread", 0);  // Don't override if already set
  setenv("OMP_PLACES", "threads", 0);
#endif
}

}  // namespace

/**
 * @brief Main entry point for the FDTD simulation program.
 *
 * Parses command-line arguments, validates configuration, initializes
 * the simulation environment (including optional Kokkos parallel
 * execution framework), and runs the FDTD simulation.
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on success, 1 on error
 */
int main(int argc, char* argv[])
{
  g_start_init_time = system_clock::now();

  // Configure parallel execution environment
  ConfigureParallelEnvironment();

#ifdef USE_KOKKOS
  // RAII: Kokkos will be automatically finalized when this object goes out of
  // scope
  KokkosScope kokkos_scope(argc, argv);
#endif

  try
  {
    // Set up command-line option parser
    cxxopts::Options options("FDTD Proxy", "FDTD acoustic wave simulation");
    options.allow_unrecognised_options();  // Allow Kokkos flags to pass through

    // Bind FDTD options to CLI parser
    fdtd::options::FdtdOptions opt;
    fdtd::options::FdtdOptions::BindCli(options, opt);

    // Parse command-line arguments
    auto result = options.parse(argc, argv);

    // Display help if requested
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      return 0;
    }

    // Validate configuration options
    opt.Validate();

    // Initialize and run FDTD simulation
    fdtd::FdtdProxy fd_sim(opt);
    Compute(fd_sim);

    // Report total execution time
    const auto total_duration =
        duration_cast<nanoseconds>(system_clock::now() - g_start_init_time);
    std::cout << "\n";
    PrintTiming("Total execution time", total_duration);

    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "\nError: " << e.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "\nFatal error: Unknown exception" << std::endl;
    return 1;
  }
  // KokkosScope destructor automatically calls Kokkos::finalize() here
}
