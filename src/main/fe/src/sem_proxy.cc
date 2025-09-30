//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <sem_solver.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

using namespace SourceAndReceiverUtils;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;

  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const float sponge_size[3] = {spongex, spongey, spongez};
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  domain_size_[0] = opt.lx;
  domain_size_[1] = opt.ly;
  domain_size_[2] = opt.lz;

  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;

  bool isModelOnNodes = opt.isModelOnNodes;

  const SolverFactory::methodType methodType = getMethod(opt.method);
  const SolverFactory::implemType implemType = getImplem(opt.implem);
  const SolverFactory::meshType meshType = getMesh(opt.mesh);

  if (meshType == SolverFactory::Struct)
  {
    int ex = nb_elements_[0];
    float lx = domain_size_[0];
    int elem_sizex = lx / ex;
    float ly = domain_size_[1];
    int ey = nb_elements_[1];
    int elem_sizey = ly / ey;
    int ez = nb_elements_[2];
    float lz = domain_size_[2];
    int elem_sizez = lz / ez;

    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, elem_sizex, ey, elem_sizey, ez, elem_sizez, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, elem_sizex, ey, elem_sizey, ez, elem_sizez, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, elem_sizex, ey, elem_sizey, ez, elem_sizez, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == SolverFactory::Unstruct)
  {
    int ex = nb_elements_[0];
    int ey = nb_elements_[1];
    int ez = nb_elements_[2];
    float lx = domain_size_[0];
    float ly = domain_size_[1];
    float lz = domain_size_[2];

    model::CartesianParams<float, int> param(order, ex, ey, ez, lx, ly, lz,
                                             isModelOnNodes);
    model::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // time parameters
  if (opt.autodt)
  {
    float cfl_factor = (order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;

  m_solver =
      SolverFactory::createSolver(methodType, implemType, meshType, order);
  m_solver->computeFEInit(*m_mesh, sponge_size, opt.surface_sponge,
                          opt.taper_delta);

  initFiniteElem();

  // snapshots settings
  is_snapshots_ = opt.snapshots;
  if (is_snapshots_)
  {
    snap_time_interval_ = opt.snap_time_interval;
    snap_folder_ = opt.snap_folder;
    std::filesystem::create_directories(snap_folder_);
  }

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;

  if (is_snapshots_)
  {
    std::cout << "Snapshots enable every " << snap_time_interval_
              << " iteration." << std::endl;
    std::cout << "Saved in " << snap_folder_ << " folder." << std::endl;
  }
}

void SEMproxy::run()
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverData solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement, rhsWeights);

  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;

    startOutputTime = system_clock::now();

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputPnValues(indexTimeSample, i1, rhsElement[0], pnGlobal);
    }

    // Save slice in dat format
    if (is_snapshots_ && indexTimeSample % snap_time_interval_ == 0)
    {
      saveSnapshot(indexTimeSample);
    }

    // Save pressure at receiver
    const int order = m_mesh->getOrder();

    float varnp1 = 0.0;
    for (int i = 0; i < order + 1; i++)
    {
      for (int j = 0; j < order + 1; j++)
      {
        for (int k = 0; k < order + 1; k++)
        {
          int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
          int globalNodeOnElement =
              i + j * (order + 1) + k * (order + 1) * (order + 1);
          varnp1 +=
              pnGlobal(nodeIdx, i2) * rhsWeightsRcv(0, globalNodeOnElement);
        }
      }
    }

    pnAtReceiver(0, indexTimeSample) = varnp1;

    swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;

    totalOutputTime += system_clock::now() - startOutputTime;
  }

  // Output receiver
  saveReceiver();

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

void SEMproxy::saveSnapshot(int timeSample) const
{
  std::stringstream filename;
  filename << snap_folder_ << "/slice" << timeSample << ".dat";
  std::string str_filename = filename.str();

#ifdef USE_KOKKOS
  auto subview = Kokkos::subview(pnGlobal, Kokkos::ALL, i1);
  int middle_z = nb_nodes_[2] / 2;
  int slice_start = middle_z * nb_nodes_[0] * nb_nodes_[1];
  int slice_end = slice_start + (nb_nodes_[0] * nb_nodes_[1]);
  auto xy_slice =
      Kokkos::subview(subview, Kokkos::make_pair(slice_start, slice_end));
  FENCE
#else
  std::vector<float> column_data = pnGlobal.getColumn(i1);

  vectorReal subview(column_data.size());
  for (int i = 0; i < column_data.size(); ++i)
  {
    subview[i] = column_data[i];
  }

  int middle_z = nb_nodes_[2] / 2;
  int slice_start = middle_z * nb_nodes_[0] * nb_nodes_[1];
  int slice_end = slice_start + (nb_nodes_[0] * nb_nodes_[1]);

  vectorReal xy_slice(slice_end - slice_start);
  for (int i = 0; i < slice_end - slice_start; ++i)
  {
    xy_slice[i] = subview[slice_start + i];
  }
#endif  // USE_KOKKOS
  saveSlice(xy_slice, nb_nodes_[0], nb_nodes_[1], str_filename);
  FENCE
}

void SEMproxy::saveReceiver() const
{
  std::stringstream filename;
  filename << "receiver.dat";
  std::string str_filename = filename.str();
  std::ofstream outfile(str_filename);
  for (int i = 0; i < num_sample_; i++)
  {
    outfile << i * dt_ << " " << pnAtReceiver(0, i) << std::endl;
  }
  outfile.close();
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");

  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(1, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      1, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
  pnAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "pnAtReceiver");
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // std::cout << "All source are currently are coded on element 50." <<
  // std::endl;
  std::cout << "All source are currently are coded on middle element."
            << std::endl;
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  int lx = domain_size_[0];
  int ly = domain_size_[1];
  int lz = domain_size_[2];

  // Get source element index

  int source_index = floor((src_coord_[0] * ex) / lx) +
                     floor((src_coord_[1] * ey) / ly) * ex +
                     floor((src_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the sourc element
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh->getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElement[0], i, j, k);
        cornerCoords[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(num_sample_, dt_, f0, sourceOrder);
  for (int j = 0; j < num_sample_; j++)
  {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  // Receiver computation

  int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElementRcv[i] = receiver_index;
  }

  // Get coordinates of the corners of the receiver element
  float cornerCoordsRcv[8][3];
  I = 0;
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
        cornerCoordsRcv[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoordsRcv[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoordsRcv[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }
}

std::string formatSnapshotFilename(int id, int width = 5)
{
  std::ostringstream oss;
  // oss << "snapshot" << std::setw(width) << std::setfill('0') << id;
  oss << "snapshot" << id;
  return oss.str();
}

/**
 * Save slice in matrix formating
 * Format: space-separated matrix with blank lines between rows for 3D plotting
 */
void SEMproxy::saveSlice(const VECTOR_REAL_VIEW& host_slice, int sizex,
                         int sizey, const std::string& filepath) const
{
  std::ofstream file(filepath);
  file << std::fixed << std::setprecision(6);
  file << sizex << "\n" << sizey << "\n";

  for (int i = 0; i < sizex * sizey; ++i)
  {
    file << host_slice[i];
    if ((i + 1) % sizex == 0)
    {
      file << "\n";  // New line every sizex elements
    }
    else
    {
      file << " ";
    }
  }
  file.close();
}

SolverFactory::implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return SolverFactory::MAKUTU;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

SolverFactory::meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return SolverFactory::Struct;
  if (meshArg == "ucartesian") return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

SolverFactory::methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return SolverFactory::SEM;
  if (methodArg == "dg") return SolverFactory::DG;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}
