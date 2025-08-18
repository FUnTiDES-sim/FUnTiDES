//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "cartesian_sem_mesh.h"
#include "solverFactory.hpp"
#include "SEMsolver.hpp"
#include <cxxopts.hpp>
#include "SEMproxy.hpp"
#include "dataType.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>


SEMproxy::SEMproxy(const SemProxyOptions& opt) {
  const int order = opt.order;

  const int ex = opt.ex;
  const int ey = opt.ey;
  const int ez = opt.ez;

  const float lx = opt.lx;
  const float ly = opt.ly;
  const float lz = opt.lz;

  const SolverFactory::methodType methodType = getMethod( opt.method );
  const SolverFactory::implemType implemType = getImplem( opt.implem );
  const SolverFactory::meshType meshType = getMesh( opt.mesh );

  using Base = mesh_base::BaseMesh<float,int>;
  using CartM1 = CartesianSEMmesh<float,int,1>;
  using CartM2 = CartesianSEMmesh<float,int,2>;
  using CartM3 = CartesianSEMmesh<float,int,3>;
  using UnstrM1 = CartesianUnstructMesh<float,int,1>;
  using UnstrM2 = CartesianUnstructMesh<float,int,2>;
  using UnstrM3 = CartesianUnstructMesh<float,int,3>;
  using MeshVar = std::variant<CartM1, CartM2, CartM3, UnstrM1, UnstrM2, UnstrM3>;
  using MeshParams = std::variant<CartesianParams<float,int>, CartesianUnstructParams<float,int>>;

  MeshParams params;
  if (meshType == SolverFactory::CARTESIAN) {
    params.emplace<CartesianParams<float,int>>(order, ex, ey, ez, lx, ly, lz);
  } else if (meshType == SolverFactory::UNSTRUCT_CARTESIAN) {
    params.emplace<CartesianUnstructParams<float,int>>(order, ex, ey, ez, lx, ly, lz);
  }

  MeshVar mesh;
  std::visit([&mesh, meshType, order](const auto& p) {
    if (meshType == SolverFactory::CARTESIAN) {
      switch (order) {
        case 1: mesh.emplace<CartM1>(p); break;
        case 2: mesh.emplace<CartM2>(p); break;
        case 3: mesh.emplace<CartM3>(p); break;
        default: throw std::invalid_argument("Order must be within [1, 3]");
      }
    } else if (meshType == SolverFactory::UNSTRUCT_CARTESIAN) {
      switch (order) {
        case 1: mesh.emplace<UnstrM1>(p); break;
        case 2: mesh.emplace<UnstrM2>(p); break;
        case 3: mesh.emplace<UnstrM3>(p); break;
        default: throw std::invalid_argument("Order must be within [1, 3]");
      }
    }
  }, params);

  if (meshType == SolverFactory::CARTESIAN) {
    CartesianParams<float,int> cartParams(order, ex, ey, ez, lx, ly, lz);
    switch (order) {
      case 1: mesh.emplace<CartM1>(cartParams); break;
      case 2: mesh.emplace<CartM2>(cartParams); break;
      case 3: mesh.emplace<CartM3>(cartParams); break;
      default: throw std::invalid_argument("Order must be within [1, 3]");
    }
  } else if (meshType == SolverFactory::UNSTRUCT_CARTESIAN) {
    CartesianUnstructParams<float,int> unstrParams(order, ex, ey, ez, lx, ly, lz);
    switch (order) {
      case 1: mesh.emplace<UnstrM1>(unstrParams); break;
      case 2: mesh.emplace<UnstrM2>(unstrParams); break;
      case 3: mesh.emplace<UnstrM3>(unstrParams); break;
      default: throw std::invalid_argument("Order must be within [1, 3]");
    }
  }
  // using Base = mesh_base::BaseMesh<float,int>;
  // using CartM1 = CartesianSEMmesh<float,int,1>;
  // using CartM2 = CartesianSEMmesh<float,int,2>;
  // using CartM3 = CartesianSEMmesh<float,int,3>;
  // using UnstrM1 = CartesianUnstructMesh<float,int,1>;
  // using UnstrM2 = CartesianUnstructMesh<float,int,2>;
  // using UnstrM3 = CartesianUnstructMesh<float,int,3>;
  // using MeshVar = std::variant<CartM1, CartM2, CartM3, UnstrM1, UnstrM2, UnstrM3>;
  // using ParamsPtr = std::unique_ptr<BaseParams>;

  // template<typename IndexType, typename ValueType>
  // using MeshParams = std::variant<
  //           CartesianParams<IndexType, ValueType>,
  //           UnstructuredParams<IndexType, ValueType>>;

  // ParamsPtr params;
  // if (meshType == MeshType::Cartesian) {
  //   params = std::make_unique<CartParams>(order, ex, ey, ez, lx, ly, lz);
  // } else {
  //   params = std::make_unique<UnstrParams>(order, nodes, connectivity);
  // }

  // MeshVar mesh;
  // CartesianParams<int,float> params{order, ex, ey, ez, lx, ly, lz};

  // if (meshType == SolverFactory::CARTESIAN) {
  //   switch (order) {
  //     case 1: mesh.emplace<CartM1>(params); break;
  //     case 2: mesh.emplace<CartM2>(params); break;
  //     case 3: mesh.emplace<CartM3>(params); break;
  //     default: throw std::invalid_argument("Order must be within [1, 3]");
  //   }
  // } else if (meshType == SolverFactory::UNSTRUCT_CARTESIAN) {
  //   switch (order) {
  //     case 1: mesh.emplace<UnstrM1>(params); break;
  //     case 2: mesh.emplace<UnstrM2>(params); break;
  //     case 3: mesh.emplace<UnstrM3>(params); break;
  //     default: throw std::invalid_argument("Order must be within [1, 3]");
  //   }
  // } else {
  //   throw std::invalid_argument("Invalid mesh type");
  // }

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType, order);

  // If solver has overloads/templates for each mesh:
  std::visit([&](auto& m){
    m_solver->computeFEInit(m);   // overload or template on T
  }, mesh);

  myMesh =
    std::visit([](auto& m) -> const Base* { return static_cast<Base*>(&m); }, mesh);


  nb_elements[0] = ex;
  nb_elements[1] = ey;
  nb_elements[2] = ez;

  nb_nodes[0] = ex * order + 1;
  nb_nodes[1] = ey * order + 1;
  nb_nodes[2] = ez * order + 1;

  initFiniteElem();

  std::cout << "Starting simulation with Cartesian Mesh of size "
            << "(" << ex << ',' << ey << ',' << ez << ')' << std::endl;
  std::cout << "Number of node is " << myMesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << myMesh->getNumberOfElements() << std::endl;
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
  if ( meshArg == "cartesian" ) return SolverFactory::CARTESIAN;
  if ( meshArg == "ucartesian" ) return SolverFactory::UNSTRUCT_CARTESIAN;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument( "Mesh type does not follow any valid type." );
}

SolverFactory::methodType SEMproxy::getMethod ( string methodArg )
{
  if ( methodArg == "sem" ) return SolverFactory::SEM;
  if ( methodArg == "dg" ) return SolverFactory::DG;

  throw std::invalid_argument( "Method type does not follow any valid type." );
}
