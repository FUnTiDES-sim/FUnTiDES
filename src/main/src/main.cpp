//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "SEMproxy.hpp"

#ifdef USE_CALIPER
#include "caliperUtils.hpp"
#include <caliper/cali-manager.h>
#endif // USE_CALIPER

#ifdef USE_EZV
#include "ezvLauncher.hpp"
#include <thread>
#endif // USE_EZV

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
  cout << "Elapsed Initial Time : "
       << (startRunTime - startInitTime).count() / 1E9 << " seconds." << endl;
  cout << "Elapsed Compute Time : "
       << (system_clock::now() - startRunTime).count() / 1E9 << " seconds."
       << endl;
};

void compute_loop(SEMproxy semsim) { compute(semsim); }

int main(int argc, char *argv[]) {

  startInitTime = system_clock::now();

#ifdef USE_EZV
  init_ezv();
#endif // USE_EZV

#ifdef USE_KOKKOS
  cout << "Using Kokkos" << endl;
  Kokkos::initialize(argc, argv);
  {
#endif

    cout << "\n+================================= " << endl;
    cout << "| Initializing SEM Application ... " << endl;
    cout << "+================================= \n" << endl;

    // SEMproxy semsim(argc, argv);
    constexpr int order = 3;
    SEMproxy semsim(100, 100, 100, 2000., order);

    semsim.initFiniteElem();

#ifdef USE_EZV
    ezv_init_mesh(semsim, &mesh);
    std::thread compute_thread(compute_loop, semsim);
    ezv_loop();

    if (compute_thread.joinable())
      compute_thread.join();
#else
  compute_loop(semsim);
#endif // USE_EZV

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExe Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
