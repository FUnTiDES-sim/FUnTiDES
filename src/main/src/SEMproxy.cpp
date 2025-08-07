//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "SEMproxy.hpp"

#include "solverFactory.hpp"
#include "SEMsolver.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

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

  nb_elements[0] = ex;
  nb_elements[1] = ey;
  nb_elements[2] = ez;

  nb_nodes[0] = ex * order + 1;
  nb_nodes[1] = ey * order + 1;
  nb_nodes[2] = ez * order + 1;

  std::cout << "Starting simulation with Cartesian Mesh of size "
            << "(" << ex << ',' << ey << ',' << ez << ')' << std::endl;
  std::cout << "Size of the domain is "
            << "(" << lx << ',' << ly << ',' << lz << ')' << std::endl;

  Mesh cartesianMesh(ex, ey, ez, lx, ly, lz, order);
  myMesh = cartesianMesh;
  m_solver->computeFEInit( cartesianMesh );

}

// Run the simulation.
void SEMproxy::run()
{

  time_point<steady_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverData solverData(  i1, i2, myRHSTerm, pnGlobal, rhsElement, rhsWeights );

  for (int indexTimeSample = 0; indexTimeSample < m_numSamples; indexTimeSample++)
  {
    startComputeTime = steady_clock::now();
    m_solver->computeOneStep( m_dt, indexTimeSample, solverData );

    totalComputeTime += steady_clock::now() - startComputeTime;

    startOutputTime = steady_clock::now();

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputPnValues(indexTimeSample, i1, rhsElement[0], pnGlobal);
    }

    if (indexTimeSample % 10 == 0)
    {
      std::stringstream filename;
      filename << "slice" << indexTimeSample << ".dat";
      std::string str_filename = filename.str();

      auto subview = Kokkos::subview(pnGlobal, Kokkos::ALL, i1);
      auto slice = myMesh.extractXYSlice(subview, nb_nodes[0], nb_nodes[0]/2);
      saveSlice(slice, nb_nodes[0], str_filename);
    }

    swap(i1, i2);
    swap( solverData.m_i1, solverData.m_i2 );

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
// void SEMproxy::getMeshInfo() {
//   // get information from mesh
//   m_numberOfNodes = myMesh.getNumberOfNodes();
//   m_numberOfElements = myMesh.getNumberOfElements();
//   m_numberOfPointsPerElement = myMesh.getNumberOfPointsPerElement();
//   m_numberOfInteriorNodes = myMesh.getNumberOfInteriorNodes();
// }

// Initialize arrays
void SEMproxy::init_arrays() {
  cout << "Allocate host memory for source and pressure values ..." << endl;
  myRHSTerm = allocateArray2D<arrayReal>(m_numberOfRHS,
                                         m_numSamples, "RHSTerm");
  rhsElement = allocateVector<vectorInt>(m_numberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(m_numberOfRHS, myMesh.getNumberOfPointsPerElement(), "RHSWeight");
  pnGlobal = allocateArray2D<arrayReal>(myMesh.getNumberOfNodes(), 2, "pnGlobal");
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(m_numberOfRHS, 3, "RHSLocation");
  // set number of rhs and location
  myRHSLocation(0, 0) = myMesh.domainSize(0) / 2;
  myRHSLocation(0, 1) = myMesh.domainSize(1) / 2;
  myRHSLocation(0, 2) = myMesh.domainSize(2) / 2;
  cout << "\nSource 1 location: " << myRHSLocation(0, 0) << ", "
       << myRHSLocation(0, 1) << ", " << myRHSLocation(0, 2) << endl;
  cout << "Corresponding to element id "
       << myMesh.elementFromCoordinate(myRHSLocation(0,0), myRHSLocation(0,1),myRHSLocation(0,2)) << endl;
  for (int i = 0; i < m_numberOfRHS; i++)
  {
    // extract element number for current rhs
    rhsElement[i] = myMesh.elementFromCoordinate(myRHSLocation(i,0), myRHSLocation(i,1),myRHSLocation(i,2));
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
  // Setting the weight for source ponderation on node.
  for (int i = 0; i < m_numberOfRHS; i++)
  {
    for (int j = 0; j < myMesh.getNumberOfPointsPerElement(); j++)
    {
      rhsWeights(i, j) = 1.0/ myMesh.getNumberOfPointsPerElement();
    }
  }
}

std::string formatSnapshotFilename(int id, int width = 5) {
  std::ostringstream oss;
  // oss << "snapshot" << std::setw(width) << std::setfill('0') << id;
  oss << "snapshot" << id;
  return oss.str();
}

/**
 * Save slice in gnuplot matrix format (default - best for gnuplot)
 * Format: space-separated matrix with blank lines between rows for 3D plotting
 */
void SEMproxy::saveSlice(const VECTOR_REAL_VIEW& host_slice,
               int size, const std::string& filepath) {
    std::ofstream file(filepath);

    file << std::fixed << std::setprecision(6);

    file << size << "\n" << size << "\n";

    // Standard matrix format - works with 'plot "file" matrix with image'
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            file << host_slice(y * size + x);
            file << " ";
        }
        // file << "\n";
    }

    file.close();
}
