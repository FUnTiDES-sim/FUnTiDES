//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "SEMproxy.hpp"
#ifdef USE_SHIVA
#include <functions/LagrangeBasis.hpp>
#include <functions/Quadrature.hpp>
#include <functions/Spacing.hpp>
#include <common/ShivaMacros.hpp>
#include <common/pmpl.hpp>
#include <common/types.hpp>
#endif


int main( int argc, char *argv[] )
{

using namespace shiva;
using namespace shiva::functions;
using BasisType = LagrangeBasis< double, 5, GaussLobattoSpacing >;
static constexpr double coord = 0.3;
constexpr double    value = BasisType::template value< 1 >( coord );
std::cout<< "simple string------------------> value="<<value<<std::endl;
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
