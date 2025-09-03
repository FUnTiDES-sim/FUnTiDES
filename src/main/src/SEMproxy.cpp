//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "solverFactory.hpp"
#include "SEMsolver.hpp"
#include <cxxopts.hpp>
#include "SEMproxy.hpp"
#include "dataType.hpp"
#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>


SEMproxy::SEMproxy(const SemProxyOptions& opt) {
  const int order = opt.order;
  const int ex = opt.ex;
  int ey = opt.ey;
  int ez = opt.ez;
  const float lx = opt.lx;
  float ly = opt.ly;
  float lz = opt.lz;

  const SolverFactory::methodType methodType = getMethod( opt.method );
  const SolverFactory::implemType implemType = getImplem( opt.implem );
  const SolverFactory::meshType meshType = getMesh( opt.mesh );

  if (meshType == SolverFactory::Struct) {
    switch(order) {
      case 1: {
        model_builder::CartesianStructBuilder<float, int, 1> builder;
        m_mesh_storage = builder.getModel(ex, lx/ex);
        break;
      }
      case 2: {
        model_builder::CartesianStructBuilder<float, int, 2> builder;
        m_mesh_storage = builder.getModel(ex, lx/ex);
        break;
      }
      case 3: {
        model_builder::CartesianStructBuilder<float, int, 3> builder;
        m_mesh_storage = builder.getModel(ex, lx/ex);
        break;
      }
      default:
        throw std::runtime_error("Order other than 1 2 3 is not supported (semproxy)");
    }
    ey = ex;
    ez = ex;
    ly = lx;
    lz = lx;
  }
  else if (meshType == SolverFactory::Unstruct) {
    model_builder::CartesianParams<float, int> param(order, ex, ey, ez,  lx, ly, lz);
    model_builder::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh_storage = builder.getModel();
  }
  else {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  m_mesh = std::visit([](auto& mesh) -> model::ModelApi<float, int>* {
      return static_cast<model::ModelApi<float, int>*>(&mesh);
    }, m_mesh_storage);

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType, order);
  m_solver->computeFEInit(*m_mesh);

  initFiniteElem();

  std::cout << "Starting simulation with Cartesian Mesh of size "
            << "(" << ex << ',' << ey << ',' << ez << ')' << std::endl;
  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements() << std::endl;
  std::cout << "Size of the domain is "
            << "(" << lx << ',' << ly << ',' << lz << ')' << std::endl;
  std::cout << "Launching the Method " << opt.method
            << ", the implementation " << opt.implem
            << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
}

void SEMproxy::run() {
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverData solverData(  i1, i2, myRHSTerm, pnGlobal, rhsElement, rhsWeights);

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

    // TODO: redo snapshot
    // if (indexTimeSample % 10 == 0)
    // {
    //   std::stringstream filename;
    //   filename << "slice" << indexTimeSample << ".dat";
    //   std::string str_filename = filename.str();

    //   auto subview = Kokkos::subview(pnGlobal, Kokkos::ALL, i1);
    //   auto slice = m_mesh->extractXYSlice(subview, nb_nodes[0], nb_nodes[0]/2);
    //   saveSlice(slice, nb_nodes[0], str_filename);
    // }

    swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;


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
  rhsWeights = allocateArray2D<arrayReal>(myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
}

// Initialize sources
void SEMproxy::init_source() {
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  std::cout << "All source are currently are coded on element 50." << std::endl;
  for (int i = 0; i < 1; i++) {
    // TODO: Make a better source init
    rhsElement[i] = 50;
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
    for (int j = 0; j < m_mesh->getNumberOfPointsPerElement(); j++)
    {
      rhsWeights(i, j) = 1.0 / m_mesh->getNumberOfPointsPerElement();
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
 * Save slice in matrix formating
 * Format: space-separated matrix with blank lines between rows for 3D plotting
 */
void SEMproxy::saveSlice(const VECTOR_REAL_VIEW& host_slice,
               int size, const std::string& filepath)
{
    std::ofstream file(filepath);

    file << std::fixed << std::setprecision(6);
    file << size << "\n" << size << "\n";

    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            file << host_slice[y * size + x];
            file << " ";
        }
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

SolverFactory::meshType SEMproxy::getMesh ( string meshArg )
{
  if ( meshArg == "cartesian" ) return SolverFactory::Struct;
  if ( meshArg == "ucartesian" ) return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument( "Mesh type does not follow any valid type." );
}

SolverFactory::methodType SEMproxy::getMethod ( string methodArg )
{
  if ( methodArg == "sem" ) return SolverFactory::SEM;
  if ( methodArg == "dg" ) return SolverFactory::DG;

  throw std::invalid_argument( "Method type does not follow any valid type." );
}
