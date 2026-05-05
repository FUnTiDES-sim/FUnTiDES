//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************
#include <cstdlib>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdlib.h>

#include "sem_proxy.h"
#include "sem_proxy_options.h"

void init_mpi(int argc, char **argv, int *rank, int *size) {
#ifdef USE_MPI
  std::cout << "Initialize MPI..." << std::endl;

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, size);
#else
  std::cout << "No MPI initialization." << std::enld;
#endif
}

void finalize_mpi() {
#ifdef USE_MPI
  std::cout << "Initialize MPI..." << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  std::cout << "No MPI involved. No de-init needed." << std::endl;
#endif
}

SemProxyOptions parse_option(int argc, char **argv) {
  cxxopts::Options options("SEM Proxy", "Runs the SEM simulation.");
  options.allow_unrecognised_options();  // lets Kokkos flags pass
  options.add_options()("h,help", "Print help message");

  SemProxyOptions opt;
  SemProxyOptions::bind_cli(options, opt);

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  try {
    opt.validate();
  } catch (const std::exception &e) {
    std::cerr << "Invalid options: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }

  return opt;
}

int main(int argc, char **argv) {
  int rank, size;
  init_mpi(argc, argv, &rank, &size);

  setenv("OMP_PROC_BIND", "spread", 1);
  setenv("OMP_PLACES", "threads", 1);
  Kokkos::initialize(argc, argv);
  {
    auto opt = parse_option(argc, argv);
    SEMproxy semsim(opt);
    std::cout << "Launching simulation." << std::endl;
    semsim.run();
    std::cout << "Ending simulation." << std::endl;
  }
  Kokkos::finalize();
  finalize_mpi();

  return (0);
}
