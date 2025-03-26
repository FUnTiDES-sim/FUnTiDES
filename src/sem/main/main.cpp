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
#endif // USE_EZV

int main(int argc, char *argv[]) {

  time_point<system_clock> startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  cout << "Using Kokkos" << endl;
  Kokkos::initialize(argc, argv);
  {
#endif

#ifdef USE_CALIPER
    cali::ConfigManager mgr;
    int caliperInitRet = launch_caliper_ctx(argc, argv, mgr);
    CALI_CXX_MARK_FUNCTION;
#endif // USE_CALIPER

    cout << "\n+================================= " << endl;
    cout << "| Initializing SEM Application ... " << endl;
    cout << "+================================= \n" << endl;

    SEMproxy semsim(argc, argv);

#ifdef USE_EZV
    static mesh3d_obj_t mesh;
    init_ezv(semsim, mesh);
#endif // USE_EZV

    semsim.initFiniteElem();

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

#ifdef USE_CALIPER
    mgr.flush();
#endif // USE_CALIPER

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExe Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
