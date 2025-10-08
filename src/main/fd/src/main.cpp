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

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include "fdtd_options.h"
#include "fdtd_proxy.h"

using std::chrono::system_clock;
using std::chrono::time_point;

// Global start time for total execution timing
time_point<system_clock> g_start_init_time;

/**
 * @brief Executes the complete FDTD simulation workflow.
 *
 * Initializes the FDTD simulation environment, runs the time-stepping
 * loop, and reports performance metrics.
 *
 * @param fd_sim Reference to the FdtdProxy simulation object
 */
void Compute(FdtdProxy& fd_sim)
{
  // Initialize FDTD simulation
  fd_sim.InitFdtd();
  std::cout << "FDTD initialization done." << std::endl;

  // Start computation timer
  time_point<system_clock> start_run_time = system_clock::now();

  // Run simulation
  std::cout << "Starting FDTD computation..." << std::endl;
  fd_sim.Run();

  // Report timing information
  auto init_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      start_run_time - g_start_init_time);
  auto compute_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      system_clock::now() - start_run_time);

  std::cout << "Elapsed initialization time: " << init_duration.count() / 1E9
            << " seconds." << std::endl;
  std::cout << "Elapsed computation time: " << compute_duration.count() / 1E9
            << " seconds." << std::endl;
}

/**
 * @brief Wrapper function for the computation loop.
 *
 * This function provides a clean interface for the main simulation
 * execution, allowing for future extensions such as parameter sweeps
 * or ensemble runs.
 *
 * @param fd_sim Reference to the FdtdProxy simulation object
 */
void ComputeLoop(FdtdProxy& fd_sim) { Compute(fd_sim); }

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

#ifdef USE_KOKKOS
  // Configure OpenMP thread binding for optimal performance
  setenv("OMP_PROC_BIND", "spread", 1);
  setenv("OMP_PLACES", "threads", 1);
  Kokkos::initialize(argc, argv);
  {
#endif

    // Set up command-line option parser
    cxxopts::Options options("FDTD Proxy", "FDTD acoustic wave simulation");
    options.allow_unrecognised_options();  // Allow Kokkos flags to pass through

    // Bind FDTD options to CLI parser
    FdtdOptions opt;
    FdtdOptions::BindCli(options, opt);

    // Parse command-line arguments
    auto result = options.parse(argc, argv);

    // Display help if requested
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      return 0;
    }

    // Validate configuration options
    try
    {
      opt.Validate();
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error: Invalid configuration - " << e.what() << std::endl;
      return 1;
    }

    // Initialize FDTD simulation object
    FdtdProxy fd_sim(opt);

    // Execute simulation
    ComputeLoop(fd_sim);

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  // Report total execution time
  auto total_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      system_clock::now() - g_start_init_time);
  std::cout << "\nTotal execution time: " << total_duration.count() / 1E9
            << " seconds." << std::endl;

  return 0;
}
