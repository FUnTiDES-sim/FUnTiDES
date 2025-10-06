//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "fdtd_options.h"
#include "fdtd_proxy.h"

time_point<system_clock> startInitTime;

void compute(fdtd_proxy &fdsim) 
{
  // initialize FD simulation
  fdsim.init_fdtd();
  cout << "FD Initialization done." << endl;
  // start timer
  time_point<system_clock> startRunTime = system_clock::now();
  // run simulation
  cout<<"Starting FD computation ... " << endl;
  fdsim.run();
  // print timing information
  cout << "Elapsed Initial Time : "
       << (startRunTime - startInitTime).count() / 1E9 << " seconds." << endl;
  cout << "Elapsed Compute Time : "
       << (system_clock::now() - startRunTime).count() / 1E9 << " seconds."
       << endl;
};

void compute_loop(fdtd_proxy & fdsim) { compute(fdsim); }

int main(int argc, char *argv[]) {

  startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  setenv("OMP_PROC_BIND", "spread", 1);
  setenv("OMP_PLACES", "threads", 1);
  Kokkos::initialize(argc, argv);
  {
#endif

  cxxopts::Options options("FD Proxy", "Runs the FD simulation.");
  options.allow_unrecognised_options();       // lets Kokkos flags pass

  fdtd_options opt;
  fdtd_options::bind_cli(options, opt);

  auto result = options.parse(argc, argv);

  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  try { opt.validate(); }
  catch (const std::exception& e) 
  {
    // your error path (no help printing here)
    std::cerr << "Invalid options: " << e.what() << "\n";
    return 1;
  }

  // initialize FD simulation object
  fdtd_proxy fdsim(opt);
  // run simulation
  compute_loop(fdsim);

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExe Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
