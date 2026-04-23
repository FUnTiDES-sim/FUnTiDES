//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "sem_proxy.h"
#include "sem_proxy_options.h"

time_point<system_clock> startInitTime;

void compute(SEMproxy &semsim) {
  cout << "\n+================================= " << endl;
  cout << "| Running SEM Application ...      " << endl;
  cout << "+================================= \n" << endl;

  // start timer
  time_point<system_clock> startRunTime = system_clock::now();
  semsim.run();

  cout << "\n+================================= " << endl;
  cout << "| SEM Application Finished.       " << endl;
  cout << "+================================= \n" << endl;

  // print timing information
  cout << "Elapsed Initial Time : " << (startRunTime - startInitTime).count() / 1E9 << " seconds." << endl;
  cout << "Elapsed Compute Time : " << (system_clock::now() - startRunTime).count() / 1E9 << " seconds." << endl;
};

void compute_loop(SEMproxy &semsim) { compute(semsim); }

// Helper to determine local rank for GPU affinity
int get_local_rank() {
#ifdef USE_MPI
  // Check standard environment variables first
  if (const char *p = std::getenv("OMPI_COMM_WORLD_LOCAL_RANK")) return std::atoi(p);
  if (const char *p = std::getenv("MV2_COMM_WORLD_LOCAL_RANK")) return std::atoi(p);
  if (const char *p = std::getenv("SLURM_LOCALID")) return std::atoi(p);

  // Fallback to MPI split type if env vars are missing
  MPI_Comm nodeComm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nodeComm);
  int local_rank;
  MPI_Comm_rank(nodeComm, &local_rank);
  MPI_Comm_free(&nodeComm);
  return local_rank;
#else
  return 0;
#endif
}

int main(int argc, char *argv[]) {
  startInitTime = system_clock::now();

#ifdef USE_MPI
  // mpi
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  setenv("OMP_PROC_BIND", "spread", 1);
  setenv("OMP_PLACES", "threads", 1);
  Kokkos::initialize(argc, argv);
  {
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
      // your error path (no help printing here)
      std::cerr << "Invalid options: " << e.what() << "\n";
      return 1;
    }

    cout << "+==================================+" << endl;
    cout << "| Initializing SEM Application ... |" << endl;
    cout << "+==================================+\n" << endl;

    SEMproxy semsim(opt);

    compute_loop(semsim);
  }
  Kokkos::finalize();

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  cout << "Elapsed TotalExe Time : " << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n" << endl;
  return (0);
}
