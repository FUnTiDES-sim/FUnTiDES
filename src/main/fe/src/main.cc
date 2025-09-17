//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "sem_proxy.h"
#include "sem_proxy_options.h"

time_point<system_clock> startInitTime;

void compute(SEMproxy &semsim)
{
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
  cout << "Elapsed Initial Time : "
       << (startRunTime - startInitTime).count() / 1E9 << " seconds." << endl;
  cout << "Elapsed Compute Time : "
       << (system_clock::now() - startRunTime).count() / 1E9 << " seconds."
       << endl;
};

void compute_loop(SEMproxy &semsim) { compute(semsim); }

int main(int argc, char *argv[])
{
  startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  setenv("OMP_PROC_BIND", "spread", 1);
  setenv("OMP_PLACES", "threads", 1);
  Kokkos::initialize(argc, argv);
  {
#endif

    cxxopts::Options options("SEM Proxy", "Runs the SEM simulation.");
    options.allow_unrecognised_options();  // lets Kokkos flags pass

    options.add_options()("h,help", "Print help message");

    SemProxyOptions opt;
    SemProxyOptions::bind_cli(options, opt);

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    try
    {
      opt.validate();
    }
    catch (const std::exception &e)
    {
      // your error path (no help printing here)
      std::cerr << "Invalid options: " << e.what() << "\n";
      return 1;
    }

    cout << "+==================================+" << endl;
    cout << "| Initializing SEM Application ... |" << endl;
    cout << "+==================================+\n" << endl;

    SEMproxy semsim(opt);

    compute_loop(semsim);

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExe Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
