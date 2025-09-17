//************************************************************************
// FD proxy application v.0.0.1
//
// main.cpp: this main file is simply a driver
//************************************************************************

#include <fd/FDTDinit.hpp>

int main(int argc, char *argv[])
{
  time_point<system_clock> startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  Kokkos::initialize(argc, argv);
  {
#endif

    cout << "\n+==================================== " << endl;
    cout << "| Initializing FDTD Application ...       " << endl;
    cout << "+==================================== \n" << endl;

    struct FDTDGRIDS myGrids;
    struct FDTDMODELS myModels;

    // imports utility and kernel modules
    FDTDInit myInit;
    FDTDKernel myKernel;
    FDTDUtils myFDTDUtils;

    // initialize geometry
    myInit.init_geometry(argc, argv, myGrids);

    // initialize coefficients
    myInit.init_coefficients(myGrids, myModels);

    // initialize source
    myInit.init_source(myModels);

    // initialize velocity and pressure models, etc
    myInit.init_models(myGrids, myModels);

    // initialize sponge boundary
    myInit.defineSpongeBoundary(myGrids, myModels.spongeArray);

    cout << "\n+================================= " << endl;
    cout << "|  Running FDTD Application ...       " << endl;
    cout << "+================================= \n" << endl;

    // start timer
    time_point<system_clock> startRunTime = system_clock::now();

    time_point<system_clock> startAddRHS, totalAddRHS, startComputeTime,
        startOutputTime, totalComputeTime, totalOutputTime;

    // main loop for wave propagation on each time step
    for (int itSample = 0; itSample < myInit.nSamples; itSample++)
    {
      // add RHS term
      startAddRHS = system_clock::now();
      myKernel.addRHS(myGrids, itSample, myModels, myInit.i2,
                      myModels.pnGlobal);
      totalAddRHS += system_clock::now() - startAddRHS;

      // compute one step
      startComputeTime = system_clock::now();
      if (myInit.usePML)
      {
        myKernel.computeOneStepPML(myGrids, myModels, myInit.i1, myInit.i2);
      }
      else
      {
        myKernel.computeOneStepSB(myGrids, myModels, myInit.i1, myInit.i2);
      }
      totalComputeTime += system_clock::now() - startComputeTime;

      // swap wavefields
      startOutputTime = system_clock::now();
      swap(myInit.i1, myInit.i2);
      // print infos and save wavefields
      myFDTDUtils.output(myGrids, myModels.pnGlobal, itSample, myInit.i2,
                         myInit.saveSnapShots);
      totalOutputTime += system_clock::now() - startOutputTime;
    }

    cout << "\n+================================= " << endl;
    cout << "|  FDTD Application Finished.       " << endl;
    cout << "+================================= \n" << endl;

    // print timing information
    float addrhstime_ms =
        time_point_cast<microseconds>(totalAddRHS).time_since_epoch().count();
    float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                              .time_since_epoch()
                              .count();
    float outputtime_ms = time_point_cast<microseconds>(totalOutputTime)
                              .time_since_epoch()
                              .count();

    cout << "Elapsed Initialization Time : "
         << (startRunTime - startInitTime).count() / 1E9 << " seconds.\n"
         << endl;
    cout << "Elapsed ComputeLoop RunTime : "
         << (system_clock::now() - startRunTime).count() / 1E9 << " seconds."
         << endl;
    cout << "------------------------------------------------ " << endl;
    cout << "---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
         << endl;
    cout << "---- Elapsed AddRHS Time : " << addrhstime_ms / 1E6 << " seconds."
         << endl;
    cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
         << endl;
    cout << "------------------------------------------------ " << endl;

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExecution Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
