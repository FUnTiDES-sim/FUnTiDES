struct SEMinfo
{
  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;

  const int myNumberOfRHS=1;
  static constexpr int myOrderNumber_a=2;
  static constexpr int myOrderNumber_b=2;
  static constexpr int myOrderNumber_c=2 * (DIMENSION==3); //Only used if DIMENSION==3
  const float myTimeStep=0.001;
  const int nPointsPerElement = (myOrderNumber_a+1) * (myOrderNumber_b+1) * (myOrderNumber_c+1);

  const float f0=10.;
  const float myTimeMax=1.;
  const int sourceOrder=1;

  int myNumSamples=myTimeMax/myTimeStep;
  int myElementSource;

  #ifdef SEM_MESHCOLOR
  int numberMaxOfElementsByColor;
  const int numberOfColors=4;
  int numberOfElementsByColor[4];
  #endif
};
