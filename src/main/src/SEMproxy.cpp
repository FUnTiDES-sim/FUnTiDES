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
  m_solver = createSolver( physicsType, methodType, order );

  m_dt = ( cmdOptionExists(argv, argv + argc, "-dt")) ? std::stof(getCmdOption(argv, argv + argc, "-dt")) : 0.001;
  m_maxTime = ( cmdOptionExists(argv, argv + argc, "-tmax")) ? std::stof(getCmdOption(argv, argv + argc, "-tmax")) : 1.5;

  m_elementSource = ( cmdOptionExists(argv, argv + argc, "-elementSource")) ? std::stoi(getCmdOption(argv, argv + argc, "-elementSource")) : 0;
  m_f0 = ( cmdOptionExists(argv, argv + argc, "-f0")) ? std::stof(getCmdOption(argv, argv + argc, "-f0")) : 10.0;
  m_sourceOrder = ( cmdOptionExists(argv, argv + argc, "-sourceOrder")) ? std::stoi(getCmdOption(argv, argv + argc, "-sourceOrder")) : 1;
  m_numberOfRHS = ( cmdOptionExists(argv, argv + argc, "-numberOfRHS")) ? std::stoi(getCmdOption(argv, argv + argc, "-numberOfRHS")) : 1;
  m_numSamples = static_cast<int>(m_maxTime / m_dt);



  SEMmesh simpleMesh{ex, ey, ez, lx, ly, lz, order ,30, false};
  myMesh = simpleMesh;

}

// SEMproxy::SEMproxy(int ex, int ey, int ez, float lx) {
//   SEMmesh simpleMesh{ex, ey, ez, lx, lx, lx, myInfo.myOrderNumber,30, false};
//   myMesh = simpleMesh;
// }

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

  m_solver->computeFEInit( myMesh );
}

// Run the simulation.
void SEMproxy::run() 
{

#ifdef USE_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  time_point<steady_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  for (int indexTimeSample = 0; indexTimeSample < m_numSamples; indexTimeSample++) 
  {
    startComputeTime = steady_clock::now();
    m_solver->computeOneStep( m_dt, indexTimeSample, i1, i2, myRHSTerm,
                              pnGlobal, rhsElement);
    totalComputeTime += steady_clock::now() - startComputeTime;

    startOutputTime = steady_clock::now();
    m_solver->outputPnValues( indexTimeSample, i1, m_elementSource,
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
  m_numberOfNodes = myMesh.getNumberOfNodes();
  m_numberOfElements = myMesh.getNumberOfElements();
  m_numberOfPointsPerElement = myMesh.getNumberOfPointsPerElement();
  m_numberOfInteriorNodes = myMesh.getNumberOfInteriorNodes();
}

// Initialize arrays
void SEMproxy::init_arrays() {
  cout << "Allocate host memory for source and pressure values ..." << endl;
  myRHSTerm = allocateArray2D<arrayReal>(m_numberOfRHS,
                                         m_numSamples, "RHSTerm");
  rhsElement = allocateVector<vectorInt>(m_numberOfRHS, "rhsElement");
  pnGlobal = allocateArray2D<arrayReal>(m_numberOfNodes, 2, "pnGlobal");
}

// Initialize sources
void SEMproxy::init_source() 
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(m_numberOfRHS, 3, "RHSLocation");
  // set number of rhs and location
  myRHSLocation(0, 0) = 1001;
  myRHSLocation(0, 1) = 1001;
  myRHSLocation(0, 2) = 1001;
  cout << "\nSource location: " << myRHSLocation(0, 0) << ", "
       << myRHSLocation(0, 1) << ", " << myRHSLocation(0, 2) << endl;
  for (int i = 0; i < m_numberOfRHS; i++) 
  {
    // extract element number for current rhs
    rhsElement[i] = myMesh.getElementNumberFromPoints( myRHSLocation(i, 0), myRHSLocation(i, 1), myRHSLocation(i, 2));
  }

  // initialize source term
  vector<float> sourceTerm = myUtils.computeSourceTerm( m_numSamples, m_dt, m_f0, m_sourceOrder);
  for (int j = 0; j < m_numSamples; j++) {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }
  // get element number of source term
  m_elementSource = rhsElement[0];
  cout << "Element number for the source location: " << m_elementSource
       << endl
       << endl;
}

arrayReal SEMproxy::getMyRHSTerm() const { return myRHSTerm; }

void SEMproxy::setMyRHSTerm(const arrayReal &value) { myRHSTerm = value; }

arrayReal SEMproxy::getPnGlobal() const { return pnGlobal; }

void SEMproxy::setPnGlobal(const arrayReal &value) { pnGlobal = value; }

vectorInt SEMproxy::getRhsElement() const { return rhsElement; }

void SEMproxy::setRhsElement(const vectorInt &value) { rhsElement = value; }
