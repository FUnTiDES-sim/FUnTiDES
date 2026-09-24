/**
 * @file main.cc
 * @brief Command-line driver: parses options, initializes MPI and Kokkos, and runs one SEMproxy simulation.
 */

#include <cstdlib>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "sem_proxy.h"
#include "sem_proxy_options.h"

/**
 * @brief Initializes MPI when built with USE_MPI; otherwise only prints a message.
 *
 * Requests MPI_THREAD_MULTIPLE because snapshot writing runs on background threads.
 * Prints a warning if that thread level is not provided.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments.
 * @param[out] rank MPI rank of the calling process (left unchanged without USE_MPI).
 * @param[out] size Number of MPI processes (left unchanged without USE_MPI).
 */
void InitMpi(int argc, char** argv, int* rank, int* size) {
#ifdef USE_MPI
  std::cout << "Initializing MPI..." << std::endl;

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (provided < MPI_THREAD_MULTIPLE) {
    std::cout << "WARNING: MPI_THREAD_MULTIPLE not supported. Async I/O may be unstable.\n";
  }

  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, size);
#else
  std::cout << "No MPI initialization." << std::endl;
#endif
}

/**
 * @brief Waits on a barrier and finalizes MPI when built with USE_MPI; otherwise only prints a message.
 */
void FinalizeMpi() {
#ifdef USE_MPI
  std::cout << "Finalizing MPI..." << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  std::cout << "No MPI involved. No de-init needed." << std::endl;
#endif
}

/**
 * @brief Parses the command line into simulation options and validates them.
 *
 * Unrecognised options are accepted so that Kokkos flags pass through.
 * Prints the help text and exits with EXIT_SUCCESS on --help; prints the error
 * and exits with EXIT_FAILURE if validation fails.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments.
 * @return Validated simulation options.
 */
SemProxyOptions ParseOptions(int argc, char** argv) {
  cxxopts::Options options("SEM Proxy", "Runs the SEM simulation.");
  options.allow_unrecognised_options();
  options.add_options()("h,help", "Print help message");

  SemProxyOptions opt;
  SemProxyOptions::bind_cli(options, opt);

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    exit(EXIT_SUCCESS);
  }

  try {
    opt.validate();
  } catch (const std::exception& e) {
    std::cerr << "Invalid options: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }

  return opt;
}

/**
 * @brief Program entry point.
 *
 * Initializes MPI, sets the OpenMP binding environment, initializes Kokkos,
 * then builds and runs a SEMproxy. The simulation object is scoped so that it is
 * destroyed before Kokkos::finalize().
 */
int main(int argc, char** argv) {
  int rank = 0;
  int size = 1;
  InitMpi(argc, argv, &rank, &size);

  setenv("OMP_PROC_BIND", "spread", 1);
  setenv("OMP_PLACES", "threads", 1);

  Kokkos::initialize(argc, argv);
  {
    auto opt = ParseOptions(argc, argv);
    SEMproxy semsim(opt);

    std::cout << "Launching simulation." << std::endl;
    semsim.Run();
    std::cout << "Ending simulation." << std::endl;
  }
  Kokkos::finalize();

  FinalizeMpi();

  return EXIT_SUCCESS;
}
