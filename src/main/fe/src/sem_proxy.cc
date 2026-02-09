//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <boundary_synchronizer.h>
#include <cartesian_partitioner.h>
#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <variant>

#ifdef USE_MPI
#include "mpi_backend.h"
#endif
#include "sem_solver.h"
#include "topology_factory.h"

using namespace SourceAndReceiverUtils;
using namespace solver::fe;
using namespace solver::fe::enums;
using namespace solver::fe::solver_factory;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  // Init MPI parameters
  int mpi_initialized = 0;
  init_mpi(&mpi_initialized);

  init_sim_params(opt);
  init_mesh_params(opt);
  init_topology();
  init_sync();
  init_time_params(opt);

  const methodType methodType = getMethod(opt.method);
  const implemType implemType = getImplem(opt.implem);
  const meshType meshType = getMesh(opt.mesh);
  const modelLocationType modelLocation = opt.isModelOnNodes
                                              ? modelLocationType::kOnNodes
                                              : modelLocationType::kOnElements;
  const physicType physicType =
      opt.isElastic ? physicType::kElastic : physicType::kAcoustic;

  m_solver = createSolver(methodType, implemType, meshType, modelLocation,
                          physicType, opt.order);
  const model::AnisotropyType anisotropyType = getAnisotropy(opt.anisotropy);

  // Build Mesh using LOCAL parameters
  if (meshType == meshType::kStruct)
  {
    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            isModelOnNodes, isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            isModelOnNodes, isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            isModelOnNodes, isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == meshType::kUnstruct)
  {
    // Pass local params to unstructured builder (handles origin internally)
    model::CartesianUnstructBuilder<float, int> builder(m_localParams);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // Init topology
  par_topology_ =
      TopologyFactory::createFromMesh(*m_mesh, dist_ctx_.rank, dist_ctx_.size,
                                      m_localParams.origin_x, m_localParams.lx);

  // Initialize Synchronizer
  std::unique_ptr<BoundarySynchronizer::Backend> backend;
  m_syncer = (dist_ctx_.rank > 1)
                 ? std::make_unique<BoundarySynchronizer>(
                       std::make_unique<DebugBackend>(dist_ctx_.rank))
                 : std::make_unique<BoundarySynchronizer>(
                       std::make_unique<SerialBackend>());

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

  m_solver = solver_factory::createSolver(methodType, implemType, meshType,
                                          modelLocation, physicType, order);

  if (isElastic)
  {
    m_solver->setAnisotropyType(anisotropyType);

    // Initialize elasticity tensors ONLY for TTI on elements
    // (ISO and VTI are computed on-the-fly, TTI on nodes also on-the-fly)
    if (anisotropyType == model::AnisotropyType::kTTI && !isModelOnNodes)
    {
      m_mesh->initElasticityTensors(anisotropyType);
    }
  }

  // Setup Sponge Parameters
  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};

  // Note: m_solver->computeFEInit is now called in run() to pass partition
  // info. We manually call init arrays here if needed, but computeFEInit does
  // it too. For consistency with old code structure, we prep arrays now.
  initFiniteElem();

  // Prepare ADIOS2 dimensions
  // Global Dimensions
  // Note: For unstructured/general meshes, this is complex.
  // For Cartesian Struct Mesh:
  // Global Size = Global_Elements * Order + 1
  // Local Size = Local_Elements * Order + 1
  // Offset = Origin_Index (calculated from origin_x)
  size_t global_nx =
      m_localParams.global_lx / (m_localParams.global_lx / opt.ex) * opt.order +
      1;  // Approximation if uniform
  // More precise:
  size_t g_ex = static_cast<size_t>(opt.ex);
  size_t g_ey = static_cast<size_t>(opt.ey);
  size_t g_ez = static_cast<size_t>(opt.ez);

  // Just 1D decomposition on X supported by partitioner currently
  size_t g_nx = g_ex * opt.order + 1;
  size_t g_ny = g_ey * opt.order + 1;
  size_t g_nz = g_ez * opt.order + 1;
  std::vector<size_t> global_dims = {g_nx, g_ny, g_nz};

  // Calculate offset based on origin
  // origin_x is physical coordinate. We need node index.
  // dx = lx / ex.
  float dx = m_localParams.global_lx / opt.ex;
  size_t elem_offset_x =
      static_cast<size_t>(std::round(m_localParams.origin_x / dx));
  size_t node_offset_x = elem_offset_x * opt.order;

  std::vector<size_t> start_offsets = {node_offset_x, 0, 0};

  // Local Size
  // Note: ADIOS2 overlap handling.
  // Nodes are shared. If we write all local nodes, boundary nodes are written
  // twice. ADIOS2 allows this, but we must ensure we are consistent. Simple
  // approach: Write full local block.
  size_t l_nx = static_cast<size_t>(nb_nodes_[0]);
  size_t l_ny = static_cast<size_t>(nb_nodes_[1]);
  size_t l_nz = static_cast<size_t>(nb_nodes_[2]);
  std::vector<size_t> local_dims = {l_nx, l_ny, l_nz};

  // Refine Global Dims vector
  global_dims = {g_nx, g_ny, g_nz};

  io_ctrl_ = std::make_shared<SemIOController>(
      global_dims, start_offsets, local_dims, static_cast<size_t>(num_sample_),
      static_cast<size_t>(1));

  // snapshots settings
  is_snapshots_ = opt.snapshots;
  if (is_snapshots_)
  {
    snap_time_interval_ = opt.snap_time_interval;
  }

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;

  if (isElastic)
  {
    std::cout << "Anisotropy type is " << opt.anisotropy << std::endl;
  }

  if (is_snapshots_)
  {
    std::cout << "Snapshots enable every " << snap_time_interval_
              << " iteration." << std::endl;
  }
}

void SEMproxy::run()
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  bool isElastic = isElastic_;

  // Sponge params from options
  const float spongex =
      0;  // Configured earlier but variable scope issue in original
  const std::array<float, 3> sponge_size = {0, 0, 0};
  const bool surface_sponge = false;
  const float taper_delta = 0.015;

  // Initialize Solver with Partition Info & Compute Local Mass
  m_solver->computeFEInit(*m_mesh, sponge_size, surface_sponge, taper_delta);

  // Synchronize Mass Matrix (Critical for DD)
  if (par_topology_.isDistributed())
  {
    m_syncer->synchronize(m_solver->getMassMatrix(), par_topology_);
  }

  auto& M = m_solver->getMassMatrix();
  // Get the global node index of the first node of the source element
  int debugNodeIdx = m_mesh->globalNodeIndex(myElementSource, 0, 0, 0);

  if (!isElastic)
  {
    WavefieldAcoustic wavefield(pnGlobalPrev, pnGlobalCurr);
    RhsAcoustic rhs(myRHSTerm, rhsElement, rhsWeights);
    SEMsolverDataAcoustic solverData(wavefield, rhs);

    for (int indexTimeSample = 0; indexTimeSample < num_sample_;
         indexTimeSample++)
    {
      startComputeTime = system_clock::now();

      // Compute Local Forces
      m_solver->computeForces(dt_, indexTimeSample, solverData);

      // Synchronize Forces
      if (par_topology_.isDistributed())
      {
        for (int c = 0; c < m_solver->getNumComponents(); ++c)
        {
          m_syncer->synchronize(m_solver->getForceVector(c), par_topology_);
        }
      }

      m_solver->updateSolution(dt_, solverData);

      totalComputeTime += system_clock::now() - startComputeTime;
      startOutputTime = system_clock::now();

      if (indexTimeSample % 50 == 0)
      {
        m_solver->outputSolutionValues(indexTimeSample, rhsElement[0],
                                       pnGlobalPrev, "pnGlobal");
      }

      // Save slice in dat format
      if (is_snapshots_ && indexTimeSample % snap_time_interval_ == 0)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "Save snapshot." << std::endl;
        saveSnapshot(indexTimeSample, pnGlobalPrev);
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
                pnGlobalCurr(nodeIdx) * rhsWeightsRcv(0, globalNodeOnElement);
          }
        }
      }

      pnAtReceiver(0, indexTimeSample) = varnp1;

      solverData.swapWavefields();

      totalOutputTime += system_clock::now() - startOutputTime;
    }

    for (int i = 0; i < pnAtReceiver.extent(0); i++)
    {
#ifdef USE_KOKKOS
      auto subview = Kokkos::subview(pnAtReceiver, i, Kokkos::ALL());
      vectorReal subset("receiver_save", num_sample_);
      Kokkos::deep_copy(subset, subview);
#else
      auto& subview = pnAtReceiver;
      vectorReal subset(subview.extent(0) * subview.extent(1));
      for (size_t k = 0; k < subview.extent(0); ++k)
      {
        for (size_t j = 0; j < subview.extent(1); ++j)
        {
          subset[k * subview.extent(1) + j] = subview(k, j);
        }
      }
#endif
      io_ctrl_->saveReceiver(subset, src_coord_);
    }
  }
  else
  {
    WavefieldElastic wavefield(uxnGlobalPrev, uxnGlobalCurr, uynGlobalPrev,
                               uynGlobalCurr, uznGlobalPrev, uznGlobalCurr);
    RhsElastic rhs(myRHSTermx, myRHSTermy, myRHSTermz, rhsElement, rhsWeights);
    SEMsolverDataElastic solverData(wavefield, rhs);

    for (int indexTimeSample = 0; indexTimeSample < num_sample_;
         indexTimeSample++)
    {
      startComputeTime = system_clock::now();

      // Compute Local Forces
      m_solver->computeForces(dt_, indexTimeSample, solverData);

      // Synchronize Forces
      // TODO: Getting it work within semproxy
      if (par_topology_.isDistributed())
      {
        for (int c = 0; c < m_solver->getNumComponents(); ++c)
        {
          m_syncer->synchronize(m_solver->getForceVector(c), par_topology_);
        }
      }

      // Update Solution
      m_solver->updateSolution(dt_, solverData);

      totalComputeTime += system_clock::now() - startComputeTime;
      startOutputTime = system_clock::now();

      if (indexTimeSample % 50 == 0)
      {
        m_solver->outputSolutionValues(indexTimeSample, rhsElement[0],
                                       uxnGlobalPrev, "uxnGlobal");
        m_solver->outputSolutionValues(indexTimeSample, rhsElement[0],
                                       uynGlobalPrev, "uynGlobal");
        m_solver->outputSolutionValues(indexTimeSample, rhsElement[0],
                                       uznGlobalPrev, "uznGlobal");
      }

      if (is_snapshots_ && indexTimeSample % snap_time_interval_ == 0)
      {
        saveSnapshot(indexTimeSample, uxnGlobalPrev);
      }

      // Save pressure at receiver
      const int order = m_mesh->getOrder();

      float varuxnp1 = 0.0;
      float varyunp1 = 0.0;
      float varuznp1 = 0.0;
      for (int i = 0; i < order + 1; i++)
      {
        for (int j = 0; j < order + 1; j++)
        {
          for (int k = 0; k < order + 1; k++)
          {
            int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
            int globalNodeOnElement =
                i + j * (order + 1) + k * (order + 1) * (order + 1);
            varuxnp1 +=
                uxnGlobalCurr(nodeIdx) * rhsWeightsRcv(0, globalNodeOnElement);
            varyunp1 +=
                uynGlobalCurr(nodeIdx) * rhsWeightsRcv(0, globalNodeOnElement);
            varuznp1 +=
                uznGlobalCurr(nodeIdx) * rhsWeightsRcv(0, globalNodeOnElement);
          }
        }
      }

      uxnAtReceiver(0, indexTimeSample) = varuxnp1;
      uynAtReceiver(0, indexTimeSample) = varyunp1;
      uznAtReceiver(0, indexTimeSample) = varuznp1;

      solverData.swapWavefields();

      totalOutputTime += system_clock::now() - startOutputTime;
    }

    for (int i = 0; i < uxnAtReceiver.extent(0); i++)
    {
#ifdef USE_KOKKOS
      auto subview = Kokkos::subview(uxnAtReceiver, i, Kokkos::ALL());
      vectorReal subset("receiver_save", num_sample_);
      Kokkos::deep_copy(subset, subview);
#else
      auto& subview = pnAtReceiver;
      vectorReal subset(subview.extent(0) * subview.extent(1));
      for (size_t k = 0; k < subview.extent(0); ++k)
      {
        for (size_t j = 0; j < subview.extent(1); ++j)
        {
          subset[k * subview.extent(1) + j] = subview(k, j);
        }
      }
#endif  // USE_KOKKOS
      io_ctrl_->saveReceiver(subset, src_coord_);
    }
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
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;
  const auto n_nodes = m_mesh->getNumberOfNodes();
  const auto n_elements = m_mesh->getNumberOfElements();
  const auto n_points_per_element = m_mesh->getNumberOfPointsPerElement();

  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(myNumberOfRHS, n_points_per_element,
                                          "RHSWeight");

  if (!isElastic_)
  {
    myRHSTerm =
        allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
    pnGlobalCurr = allocateVector<vectorReal>(n_nodes, "pnGlobalCurr");
    pnGlobalPrev = allocateVector<vectorReal>(n_nodes, "pnGlobalPrev");
    pnAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "pnAtReceiver");
  }
  else
  {
    myRHSTermx =
        allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTermx");
    myRHSTermy =
        allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTermy");
    myRHSTermz =
        allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTermz");
    uxnGlobalCurr = allocateVector<vectorReal>(n_nodes, "uxnGlobalCurr");
    uynGlobalCurr = allocateVector<vectorReal>(n_nodes, "uynGlobalCurr");
    uznGlobalCurr = allocateVector<vectorReal>(n_nodes, "uznGlobalCurr");
    uxnGlobalPrev = allocateVector<vectorReal>(n_nodes, "uxnGlobalPrev");
    uynGlobalPrev = allocateVector<vectorReal>(n_nodes, "uynGlobalPrev");
    uznGlobalPrev = allocateVector<vectorReal>(n_nodes, "uznGlobalPrev");
    uxnAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "uxnAtReceiver");
    uynAtReceiver =
        allocateArray2D<arrayReal>(1, num_sample_, "uynAtReceiver ");
    uznAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "uznAtReceiver");
  }
  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(1, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      1, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
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

  // NOTE: In DD, we need to adjust source coordinate relative to local origin
  // to find correct local element.
  // However, since we are using local params for ex/lx calculation here,
  // we just need to ensure src_coord_ is treated correctly.
  // The current logic calculates index relative to the local mesh.
  // We need to shift the coordinate by origin_x.

  float relative_src_x = src_coord_[0] - m_localParams.origin_x;
  float relative_src_y = src_coord_[1] - m_localParams.origin_y;
  float relative_src_z = src_coord_[2] - m_localParams.origin_z;

  // Simple check: is source on this rank?
  bool sourceOnThisRank = (relative_src_x >= 0 && relative_src_x < lx);
  // Note: Y and Z not partitioned in 1D case, so check logic is simpler.

  int source_index = 0;
  if (sourceOnThisRank)
  {
    source_index = floor((relative_src_x * ex) / lx) +
                   floor((relative_src_y * ey) / ly) * ex +
                   floor((relative_src_z * ez) / lz) * ey * ex;
  }
  else
  {
    // Point to a safe dummy element (e.g. 0) but we will zero out weights?
    // Or we can rely on weight computation finding it outside.
    // For proxy simplicity, if not local, we just calculate 'an' index.
    // This logic should ideally be robust.
    source_index = 0;  // Placeholder
  }

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the source element
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
  if (!isElastic_)
  {
    for (int j = 0; j < num_sample_; j++)
    {
      myRHSTerm(0, j) = sourceTerm[j];
      if (j % 100 == 0)
        cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
    }
  }
  else
  {
    for (int j = 0; j < num_sample_; j++)
    {
      myRHSTermx(0, j) = sourceTerm[j];
      myRHSTermy(0, j) = sourceTerm[j];
      myRHSTermz(0, j) = sourceTerm[j];
      if (j % 100 == 0)
        cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
    }
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  // Compute Weights: If source is not local, weights will be ~0 or we should
  // explicit zero them
  if (sourceOnThisRank)
  {
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
  }
  else
  {
    // Zero weights if source is not on this rank
    for (int k = 0; k < m_mesh->getNumberOfPointsPerElement(); ++k)
      rhsWeights(0, k) = 0.0f;
  }

  // Receiver computation
  // Similar logic for receiver
  float relative_rcv_x = rcv_coord_[0] - m_localParams.origin_x;
  int receiver_index = floor((relative_rcv_x * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  // Clamp index to avoid segfaults during testing if receiver is out of bounds
  if (receiver_index < 0) receiver_index = 0;
  if (receiver_index >= m_mesh->getNumberOfElements()) receiver_index = 0;

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

  std::cout << "\n--- DEBUG INFO ---" << std::endl;
  std::cout << "Source Element: " << rhsElement[0] << std::endl;
  std::cout << "Source Coord: " << src_coord_[0] << " " << src_coord_[1] << " "
            << src_coord_[2] << std::endl;

  // Print Corner Coordinates of the source element
  std::cout << "Corner Coords (Node 0): " << cornerCoords[0][0] << ", "
            << cornerCoords[0][1] << ", " << cornerCoords[0][2] << std::endl;
  std::cout << "Corner Coords (Node 7): " << cornerCoords[7][0] << ", "
            << cornerCoords[7][1] << ", " << cornerCoords[7][2] << std::endl;

  // Print Calculated Weights
  std::cout << "RHS Weights: ";
  for (int k = 0; k < m_mesh->getNumberOfPointsPerElement(); ++k)
  {
    std::cout << rhsWeights(0, k) << " ";
  }
  std::cout << std::endl;
  std::cout << "------------------\n" << std::endl;
}

void SEMproxy::saveSnapshot(int timestep, VECTOR_REAL_VIEW data) const
{
#ifdef USE_KOKKOS
  auto nb_nodes = data.extent(0);

  vectorReal subset("snapshot_cpy", nb_nodes);
  // Use a parallel copy to handle the strided layout
  Kokkos::parallel_for(
      "copy_column", nb_nodes, KOKKOS_LAMBDA(int i) { subset(i) = data(i); });
  Kokkos::fence();
#else
  auto& subset = data;
#endif  // USE_KOKKOS

  io_ctrl_->saveSnapshot(subset, timestep);
}

implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return implemType::kMakutu;
  if (implemArg == "shiva") return implemType::kShiva;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return meshType::kStruct;
  if (meshArg == "ucartesian") return meshType::kUnstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return methodType::kSem;
  if (methodArg == "dg") return methodType::kDg;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

model::AnisotropyType SEMproxy::getAnisotropy(string anisotropyArg)
{
  if (anisotropyArg == "iso") return model::AnisotropyType::kIso;
  if (anisotropyArg == "vti") return model::AnisotropyType::kVTI;
  if (anisotropyArg == "tti") return model::AnisotropyType::kTTI;

  throw std::invalid_argument(
      "Anisotropy type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}

void SEMproxy::init_mpi(int* mpi_init)
{
#ifdef USE_MPI
  MPI_Initialized(mpi_init);
  if (mpi_init)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &dist_ctx_.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &dist_ctx_.size);
  }
  else
  {
    dist_ctx_.rank = 0;
    dist_ctx_.size = 1;
  }
  std::cout << "[rank " << dist_ctx_.rank << "] size " << dist_ctx_.size
            << std::endl;
#else
  dist_ctx_.rank = 0;
  dist_ctx_.size = 1;
#endif
}

void SEMproxy::init_sim_params(const SemProxyOptions& opt)
{
  // Partition Logic
  // Create Global Params
  model::CartesianParams<float, int> globalParams(
      opt.order, opt.ex, opt.ey, opt.ez, opt.lx, opt.ly, opt.lz,
      opt.isModelOnNodes, opt.isElastic);
  globalParams.origin_x = 0;  // Global start

  // Partition domain
  model::CartesianXPartitioner<float, int> partitioner;
  m_localParams =
      partitioner.partition(globalParams, dist_ctx_.rank, dist_ctx_.size);

  // Update members with LOCAL parameters for array allocation
  nb_elements_[0] = m_localParams.ex;
  nb_elements_[1] = m_localParams.ey;
  nb_elements_[2] = m_localParams.ez;
  nb_nodes_[0] = m_localParams.ex * opt.order + 1;
  nb_nodes_[1] = m_localParams.ey * opt.order + 1;
  nb_nodes_[2] = m_localParams.ez * opt.order + 1;

  // Use local dimensions for domain size check logic
  domain_size_[0] = m_localParams.lx;
  domain_size_[1] = m_localParams.ly;
  domain_size_[2] = m_localParams.lz;

  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;

  isElastic_ = opt.isElastic;

  std::cout << "Debug Print :" << std::endl;
  std::cout << "    Rank " << dist_ctx_.rank << "/" << dist_ctx_.size
            << std::endl;
  std::cout << "    Local lx " << m_localParams.lx << std::endl;
  std::cout << "    Local ly " << m_localParams.ly << std::endl;
  std::cout << "    Local lz " << m_localParams.lz << std::endl;
  std::cout << "    Local ex " << m_localParams.ex << std::endl;
  std::cout << "    Local ey " << m_localParams.ey << std::endl;
  std::cout << "    Local ez " << m_localParams.ez << std::endl;
}

void SEMproxy::init_mesh_params(const SemProxyOptions& opt)
{
  const methodType methodType = getMethod(opt.method);
  const implemType implemType = getImplem(opt.implem);
  const meshType meshType = getMesh(opt.mesh);

  // Build Mesh using LOCAL parameters
  if (meshType == meshType::kStruct)
  {
    switch (opt.order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            opt.isModelOnNodes, opt.isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            opt.isModelOnNodes, opt.isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            opt.isModelOnNodes, opt.isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == meshType::kUnstruct)
  {
    // Pass local params to unstructured builder (handles origin internally)
    model::CartesianUnstructBuilder<float, int> builder(m_localParams);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }
}

void SEMproxy::init_topology()
{
  par_topology_ =
      TopologyFactory::createFromMesh(*m_mesh, dist_ctx_.rank, dist_ctx_.size,
                                      m_localParams.origin_x, m_localParams.lx);
}

void SEMproxy::init_sync()
{
#if USE_MPI
  if (dist_ctx_.size > 1)
  {
    m_syncer = std::make_unique<BoundarySynchronizer>(
        std::make_unique<solver::fe::MPIBackend>());

    if (dist_ctx_.rank == 0)
    {
      std::cout << "MPI Enabled: Using MPIBackend for " << dist_ctx_.size
                << " ranks." << std::endl;
    }
  }
  else
  {
    // Fallback for serial
    m_syncer = std::make_unique<BoundarySynchronizer>(
        std::make_unique<SerialBackend>());
  }
#else
  m_syncer =
      std::make_unique<BoundarySynchronizer>(std::make_unique<SerialBackend>());
#endif
}

void SEMproxy::init_time_params(const SemProxyOptions& opt)
{
  if (opt.autodt)
  {
    float cfl_factor = (opt.order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;
}

void SEMproxy::display_init_msg(const SemProxyOptions& opt)
{
  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (opt.isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (opt.isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << opt.order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;

  if (is_snapshots_)
  {
    std::cout << "Snapshots enable every " << snap_time_interval_
              << " iteration." << std::endl;
  }
}
