// Driver file for the SEM proxy simulation.

#include <cstdlib>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "sem_proxy.h"
#include "sem_proxy_options.h"

/**
 * @brief Initializes the MPI environment if USE_MPI is defined.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param rank Pointer to store the MPI rank of the calling process.
 * @param size Pointer to store the total number of MPI processes.
 */
void InitMpi(int argc, char** argv, int* rank, int* size) {
#ifdef USE_MPI
  std::cout << "Initializing MPI..." << std::endl;

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, size);
#else
  std::cout << "No MPI initialization." << std::endl;
#endif
}

/**
 * @brief Finalizes the MPI environment if USE_MPI is defined.
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
 * @brief Parses command-line arguments to build simulation options.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return Parsed and validated simulation configuration options.
 */
SemProxyOptions ParseOptions(int argc, char** argv) {
  cxxopts::Options options("SEM Proxy", "Runs the SEM simulation.");
  options.allow_unrecognised_options();  // Allows Kokkos flags to pass through
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

int main(int argc, char** argv) {
  int rank = 0;
  int size = 1;
  InitMpi(argc, argv, &rank, &size);

  // Configure OpenMP thread bindings for performance
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
