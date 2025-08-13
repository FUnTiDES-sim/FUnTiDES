//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "cartesianSEMmesh.hpp"
#include "solverFactory.hpp"
#include "SEMsolver.hpp"
#include <cxxopts.hpp>
#include "SEMproxy.hpp"
#include "dataType.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>

SEMproxy::SEMproxy(int argc, char *argv[]) {
  cxxopts::Options options("SEM Proxy", "Runs the SEM simulation.");
  options.add_options()
    ("o,order", "Order of approximation.",
         cxxopts::value<int>()->default_value("2"))

    ("ex", "Number of element on X axis. For Cartesian Mesh",
         cxxopts::value<int>()->default_value("50"))
    ("ey", "Number of element on Y axis. For Cartesian Mesh",
         cxxopts::value<int>()->default_value("50"))
    ("ez", "Number of element on Z axis. For Cartesian Mesh",
         cxxopts::value<int>()->default_value("50"))

    ("lx", "Size of the Domain. X axis. For Cartesian Mesh",
         cxxopts::value<float>()->default_value("2000."))
    ("ly", "Size of the Domain. Y axis. For Cartesian Mesh",
         cxxopts::value<float>()->default_value("2000."))
    ("lz", "Size of the Domain. Z axis. For Cartesian Mesh",
         cxxopts::value<float>()->default_value("2000."))

    ("implem", "Implentation type (classic, optim, geos, shiva)",
         cxxopts::value<string>()->default_value("optim"))
    ("method", "Method type (sem, dg)",
         cxxopts::value<string>()->default_value("sem"))
    ;

  auto result = options.parse(argc, argv);

  int order = result["order"].as<int>();

  int ex = result["ex"].as<int>();
  int ey = result["ey"].as<int>();
  int ez = result["ez"].as<int>();

  float lx = result["lx"].as<float>();
  float ly = result["ly"].as<float>();
  float lz = result["lz"].as<float>();

  SolverFactory::methodType methodType = getMethod( result["method"].as<string>() );
  SolverFactory::implemType implemType = getImplem( result["implem"].as<string>() );

  m_solver = SolverFactory::createSolver( methodType, implemType, order );

  CartesianParams<int, float> params{ order, ex, ey, ez, lx, ly, lz};
  CartesianSEMmesh<float, int, 2> cartesianMesh(params);
  myMesh = &cartesianMesh;
  m_solver->computeFEInit(*myMesh);

  nb_elements[0] = ex;
  nb_elements[1] = ey;
  nb_elements[2] = ez;

  nb_nodes[0] = ex * order + 1;
  nb_nodes[1] = ey * order + 1;
  nb_nodes[2] = ez * order + 1;

  initFiniteElem();

  std::cout << "Starting simulation with Cartesian Mesh of size "
            << "(" << ex << ',' << ey << ',' << ez << ')' << std::endl;
  std::cout << "Size of the domain is "
            << "(" << lx << ',' << ly << ',' << lz << ')' << std::endl;
  std::cout << "Launching the Method " << result["method"].as<string>()
            << " and the implementation " << result["implem"].as<string>() << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
}

void SEMproxy::run() {
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverData solverData(  i1, i2, myRHSTerm, pnGlobal, rhsElement, rhsWeights );

  for (int indexTimeSample = 0; indexTimeSample < myNumSamples;
       indexTimeSample++) {
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(myTimeStep, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;

    startOutputTime = system_clock::now();

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputPnValues(indexTimeSample, i1, rhsElement[0], pnGlobal);
    }

    // if (indexTimeSample % 10 == 0)
    // {
    //   std::stringstream filename;
    //   filename << "slice" << indexTimeSample << ".dat";
    //   std::string str_filename = filename.str();

    //   auto subview = Kokkos::subview(pnGlobal, Kokkos::ALL, i1);
    //   auto slice = myMesh->extractXYSlice(subview, nb_nodes[0], nb_nodes[0]/2);
    //   saveSlice(slice, nb_nodes[0], str_filename);
    // }

    swap(i1, i2);
    totalOutputTime += system_clock::now() - startOutputTime;
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
void SEMproxy::init_arrays() {
  cout << "Allocate host memory for source and pressure values ..." << endl;
  myRHSTerm =
      allocateArray2D<arrayReal>(myNumberOfRHS, myNumSamples, "RHSTerm");
  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(myNumberOfRHS, myMesh->getNumberOfPointsPerElement(), "RHSWeight");
  pnGlobal =
      allocateArray2D<arrayReal>(myMesh->getNumberOfNodes(), 2, "pnGlobal");
}

// Initialize sources
void SEMproxy::init_source() {
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // set number of rhs and location
  myRHSLocation(0, 0) = myMesh->domainSize(0) / 2;
  myRHSLocation(0, 1) = myMesh->domainSize(1) / 2;
  myRHSLocation(0, 2) = myMesh->domainSize(2) / 2;
  cout << "\nSource 1 location: " << myRHSLocation(0, 0) << ", "
       << myRHSLocation(0, 1) << ", " << myRHSLocation(0, 2) << endl;
  cout << "Corresponding to element id "
       << myMesh->elementFromCoordinate(myRHSLocation(0,0), myRHSLocation(0,1),myRHSLocation(0,2)) << endl;
  for (int i = 0; i < 1; i++) {
    // extract element number for current rhs
    rhsElement[i] = myMesh->elementFromCoordinate(myRHSLocation(i,0), myRHSLocation(i,1),myRHSLocation(i,2));
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(myNumSamples, myTimeStep, f0, sourceOrder);
  for (int j = 0; j < myNumSamples; j++) {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }
  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;
  // Setting the weight for source ponderation on node.
  for (int i = 0; i < myNumberOfRHS; i++)
  {
    for (int j = 0; j < myMesh->getNumberOfPointsPerElement(); j++)
    {
      rhsWeights(i, j) = 1.0 / myMesh->getNumberOfPointsPerElement();
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

SolverFactory::implemType SEMproxy::getImplem ( string implemArg )
{
  if (implemArg == "classic") return SolverFactory::CLASSIC;
  if (implemArg == "optim") return SolverFactory::OPTIM;
  if (implemArg == "geos") return SolverFactory::GEOS;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument( "Implentation type does not follow any valid type." );
}

SolverFactory::methodType SEMproxy::getMethod ( string methodArg )
{
  if ( methodArg == "sem" ) return SolverFactory::SEM;
  if ( methodArg == "dg" ) return SolverFactory::DG;

  throw std::invalid_argument( "Method type does not follow any valid type." );
}
