#ifndef SEMINFO_HPP_
#define SEMINFO_HPP_

#include "SEMmacros.hpp"

#include <cmath>

  constexpr inline int intPow(int base, unsigned int exp) 
  {
    return (exp == 0) ? 1 : base * intPow(base, exp - 1);
  }

struct SEMinfo
{

  

  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;

  const int myNumberOfRHS=1;
  static constexpr int myOrderNumber=2;
  const float myTimeStep=0.001;
  static constexpr int dims=DIMENSION;
  static constexpr int nPointsPerElement = intPow((myOrderNumber+1), dims );

  const float f0=10.;
  const float myTimeMax=0.01;
  const int sourceOrder=1;

  int myNumSamples=myTimeMax/myTimeStep;
  int myElementSource;

};
#endif
