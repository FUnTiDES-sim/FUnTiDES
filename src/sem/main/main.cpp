//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "SEMproxy.hpp"

int main( int argc, char *argv[] )
{

  time_point< system_clock > startInitTime = system_clock::now();

  #ifdef USE_KOKKOS
  Kokkos::initialize( argc, argv );
  {
  #endif

  cout << "\n+================================= "<< endl;
  cout << "| Initializing SEM Application ... "<< endl;
  cout << "+================================= \n"<< endl;

  SEMproxy semsim( argc, argv );

  semsim.initFiniteElem();

  cout << "\n+================================= "<< endl;
  cout << "| Running SEM Application ...      "<< endl;
  cout << "+================================= \n"<< endl;

  // start timer
  time_point< system_clock > startRunTime = system_clock::now();
  semsim.run();

  cout << "\n+================================= "<< endl;
  cout << "| SEM Application Finished.       "<< endl;
  cout << "+================================= \n"<< endl;

  // print timing information
  cout << "Elapsed Initial Time : "<<( startRunTime - startInitTime ).count()/1E9 <<" seconds."<< endl;
  cout << "Elapsed Compute Time : "<<( system_clock::now()-startRunTime ).count()/1E9 <<" seconds."<< endl;

  #ifdef USE_KOKKOS
  }
  Kokkos::finalize();
  #endif

  cout << "Elapsed TotalExe Time : "<<( system_clock::now()-startInitTime).count()/1E9 <<" seconds.\n"<< endl;
  return (0);
}
