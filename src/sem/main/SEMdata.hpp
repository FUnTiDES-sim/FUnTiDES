#ifndef SEMINFO_HPP_
#define SEMINFO_HPP_
struct SEMinfo
{
  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;

  const int myNumberOfRHS=1;
  static constexpr int myOrderNumber=3;
  const float myTimeStep=0.001;
  const int nPointsPerElement = pow((myOrderNumber+1), DIMENSION );

  const float f0=10.;
  const float myTimeMax=0.5;
  const int sourceOrder=1;

  int myNumSamples=myTimeMax/myTimeStep;
  int myElementSource;

};
#endif
