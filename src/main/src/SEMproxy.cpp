//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "SEMproxy.hpp"

#include "solverFactory.hpp"
#ifdef USE_EZV
#include "ezvLauncher.hpp"
#endif // USE_EZV

SEMproxy::SEMproxy(int argc, char *argv[]) 
{
  int ex = ( cmdOptionExists(argv, argv + argc, "-ex")) ? std::stoi(getCmdOption(argv, argv + argc, "-ex")) : 100;
  int ey = ( cmdOptionExists(argv, argv + argc, "-ey")) ? std::stoi(getCmdOption(argv, argv + argc, "-ey")) : ex;
  int ez = ( cmdOptionExists(argv, argv + argc, "-ez")) ? std::stoi(getCmdOption(argv, argv + argc, "-ez")) : ex;

  float lx = ( cmdOptionExists(argv, argv + argc, "-lx")) ? std::stof(getCmdOption(argv, argv + argc, "-lx")) : 2000;
  float ly = ( cmdOptionExists(argv, argv + argc, "-ly")) ? std::stof(getCmdOption(argv, argv + argc, "-ly")) : lx;
  float lz = ( cmdOptionExists(argv, argv + argc, "-lz")) ? std::stof(getCmdOption(argv, argv + argc, "-lz")) : lx;

  int const physicsType = ( cmdOptionExists(argv, argv + argc, "-physics")) ? std::stoi(getCmdOption(argv, argv + argc, "-physics")) : 0;
  int const methodType = ( cmdOptionExists(argv, argv + argc, "-method")) ? std::stoi(getCmdOption(argv, argv + argc, "-method")) : 1;
  int const order = ( cmdOptionExists(argv, argv + argc, "-order")) ? std::stoi(getCmdOption(argv, argv + argc, "-order")) : 2;
  mySolver = createSolver( physicsType, methodType, order );



  SEMmesh simpleMesh{ex, ey, ez, lx, ly, lz, myInfo.myOrderNumber,30, false};
  myMesh = simpleMesh;

}

SEMproxy::SEMproxy(int ex, int ey, int ez, float lx) {
  SEMmesh simpleMesh{ex, ey, ez, lx, lx, lx, myInfo.myOrderNumber,30, false};
  myMesh = simpleMesh;
}

// Initialize the simulation.
void SEMproxy::initFiniteElem() {
#ifdef USE_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  // get information from mesh
  getMeshInfo();

  // allocate arrays and vectors
  init_arrays();

  // initialize source and RHS
  init_source();

  mySolver->computeFEInit(myInfo, myMesh);
}

// Run the simulation.
void SEMproxy::run() {

#ifdef USE_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  time_point<steady_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  for (int indexTimeSample = 0; indexTimeSample < myInfo.myNumSamples; indexTimeSample++) 
  {
    startComputeTime = steady_clock::now();
    mySolver->computeOneStep(indexTimeSample,
                            myInfo.nPointsPerElement, i1, i2, myInfo, myRHSTerm,
                            pnGlobal, rhsElement);
    totalComputeTime += steady_clock::now() - startComputeTime;

    startOutputTime = steady_clock::now();
    mySolver->outputPnValues(myMesh, indexTimeSample, i1, myInfo.myElementSource,
                            pnGlobal);

    swap(i1, i2);
    totalOutputTime += steady_clock::now() - startOutputTime;
  }

  float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                            .time_since_epoch()
                            .count();
  float outputtime_ms =
      time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();

  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
       << endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
       << endl;
  cout << "------------------------------------------------ " << endl;
}

// Initialize arrays
void SEMproxy::getMeshInfo() {
  // get information from mesh
  myInfo.numberOfNodes = myMesh.getNumberOfNodes();
  myInfo.numberOfElements = myMesh.getNumberOfElements();
  myInfo.numberOfPointsPerElement = myMesh.getNumberOfPointsPerElement();
  myInfo.numberOfInteriorNodes = myMesh.getNumberOfInteriorNodes();
  myInfo.numberOfPointsPerElement = myMesh.getNumberOfPointsPerElement();
}

// Initialize arrays
void SEMproxy::init_arrays() {
  cout << "Allocate host memory for source and pressure values ..." << endl;
  myRHSTerm = allocateArray2D<arrayReal>(myInfo.myNumberOfRHS,
                                         myInfo.myNumSamples, "RHSTerm");
  rhsElement = allocateVector<vectorInt>(myInfo.myNumberOfRHS, "rhsElement");
  pnGlobal = allocateArray2D<arrayReal>(myInfo.numberOfNodes, 2, "pnGlobal");
}

// Initialize sources
void SEMproxy::init_source() 
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(myInfo.myNumberOfRHS, 3, "RHSLocation");
  // set number of rhs and location
  myRHSLocation(0, 0) = 1001;
  myRHSLocation(0, 1) = 1001;
  myRHSLocation(0, 2) = 1001;
  cout << "\nSource location: " << myRHSLocation(0, 0) << ", "
       << myRHSLocation(0, 1) << ", " << myRHSLocation(0, 2) << endl;
  for (int i = 0; i < myInfo.myNumberOfRHS; i++) 
  {
    // extract element number for current rhs
    rhsElement[i] = myMesh.getElementNumberFromPoints( myRHSLocation(i, 0), myRHSLocation(i, 1), myRHSLocation(i, 2));
  }

  // initialize source term
  vector<float> sourceTerm = myUtils.computeSourceTerm( myInfo.myNumSamples, myInfo.myTimeStep, myInfo.f0, myInfo.sourceOrder);
  for (int j = 0; j < myInfo.myNumSamples; j++) {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }
  // get element number of source term
  myInfo.myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myInfo.myElementSource
       << endl
       << endl;
}

arrayReal SEMproxy::getMyRHSTerm() const { return myRHSTerm; }

void SEMproxy::setMyRHSTerm(const arrayReal &value) { myRHSTerm = value; }

arrayReal SEMproxy::getPnGlobal() const { return pnGlobal; }

void SEMproxy::setPnGlobal(const arrayReal &value) { pnGlobal = value; }

vectorInt SEMproxy::getRhsElement() const { return rhsElement; }

void SEMproxy::setRhsElement(const vectorInt &value) { rhsElement = value; }
