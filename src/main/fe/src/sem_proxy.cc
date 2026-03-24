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

#include <cmath>
#include <cxxopts.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <variant>

#include "sem_solver.h"
#include "topology_factory.h"

using namespace SourceAndReceiverUtils;
using namespace solver::fe;
using namespace utils::enums;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;

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
  nb_nodes_[0] = m_localParams.ex * order + 1;
  nb_nodes_[1] = m_localParams.ey * order + 1;
  nb_nodes_[2] = m_localParams.ez * order + 1;

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

  bool isModelOnNodes = opt.isModelOnNodes;
  isElastic_ = opt.isElastic;
  freeSurface_ = opt.free_surface;
  cout << boolalpha;
  bool isElastic = isElastic_;

  const methodType methodType = getMethod(opt.method);
  const implemType implemType = getImplem(opt.implem);
  const meshType meshType = getMesh(opt.mesh);
  const modelLocationType modelLocation = isModelOnNodes
                                              ? modelLocationType::kOnNodes
                                              : modelLocationType::kOnElements;
  const physicType physicType =
      isElastic ? physicType::kElastic : physicType::kAcoustic;

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
            m_localParams.origin_y, m_localParams.origin_z,
            m_localParams.global_lx, m_localParams.global_ly,
            m_localParams.global_lz, m_localParams.global_origin_x,
            m_localParams.global_origin_y, m_localParams.global_origin_z);
        m_mesh = builder.getModel(freeSurface_);
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            isModelOnNodes, isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z,
            m_localParams.global_lx, m_localParams.global_ly,
            m_localParams.global_lz, m_localParams.global_origin_x,
            m_localParams.global_origin_y, m_localParams.global_origin_z);
        m_mesh = builder.getModel(freeSurface_);
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            m_localParams.ex, m_localParams.lx, m_localParams.ey,
            m_localParams.ly, m_localParams.ez, m_localParams.lz,
            isModelOnNodes, isElastic, m_localParams.origin_x,
            m_localParams.origin_y, m_localParams.origin_z,
            m_localParams.global_lx, m_localParams.global_ly,
            m_localParams.global_lz, m_localParams.global_origin_x,
            m_localParams.global_origin_y, m_localParams.global_origin_z);
        m_mesh = builder.getModel(freeSurface_);
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
    m_mesh = builder.getModel(freeSurface_);
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

  // Setup Sponge Parameters (store as members for use in run())
  sponge_size_ = {opt.boundaries_size, opt.boundaries_size,
                  opt.boundaries_size};
  surface_sponge_ = opt.surface_sponge;
  taper_delta_ = opt.taper_delta;

  // Set quality factors for attenuation (must be set before setSLSAttenuation
  // so that auto-fill of anelasticity coefficients from Q works)
  if (opt.qp > 0 || opt.qs > 0)
  {
    float qp = (opt.qp > 0) ? opt.qp : 1e9f;
    float qs = (opt.qs > 0) ? opt.qs : 1e9f;
    m_mesh->setQualityFactors(qp, qs);
    std::cout << "Quality factors set: Qp=" << qp << ", Qs=" << qs << std::endl;
  }

  // Enable SLS attenuation mechanism
  auto toView = [](const std::vector<float>& v, const char* name) {
    auto view = allocateVector<vectorReal>(v.size(), name);
    for (size_t i = 0; i < v.size(); ++i) view[i] = v[i];
    return view;
  };

  if (!opt.sls_reference_angular_frequencies.empty())
  {
    auto freqView = toView(opt.sls_reference_angular_frequencies, "slsFreqs");
    auto coeffView = toView(opt.sls_anelasticity_coefficients, "slsCoeffs");
    m_solver->setSLSAttenuation(freqView, coeffView);
  }
  else if (opt.qp > 0 || opt.qs > 0)
  {
    // Auto-enable SLS with default reference frequency based on source f0
    float omega0 = 2.0f * M_PI * f0;
    std::cout << "Auto-enabling SLS attenuation at omega0=" << omega0
              << " rad/s (f0=" << f0 << " Hz)" << std::endl;
    auto freqView = allocateVector<vectorReal>(1, "slsFreqAuto");
    freqView[0] = omega0;
    m_solver->setSLSAttenuation(freqView);
  }

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

  // Note: m_solver->computeFEInit is now called in run() to pass partition
  // info. We manually call init arrays here if needed, but computeFEInit does
  // it too. For consistency with old code structure, we prep arrays now.
  initFiniteElem();

  io_ctrl_ = std::make_shared<SemIOController>(
      static_cast<size_t>(m_mesh->getNumberOfNodes()),
      static_cast<size_t>(num_sample_), static_cast<size_t>(1));

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

  // DAS receiver setup
  if (opt.das_type == "dipole")
  {
    dasType_ = SourceAndReceiverUtils::DASType::kDipole;
    dasNumSamples_ = 2;  // dipole uses exactly 2 endpoints
  }
  else if (opt.das_type == "strain")
  {
    dasType_ = SourceAndReceiverUtils::DASType::kStrainIntegration;
    dasNumSamples_ = opt.das_samples;
  }
  else
  {
    dasType_ = SourceAndReceiverUtils::DASType::kNone;
  }

  if (dasType_ != SourceAndReceiverUtils::DASType::kNone)
  {
    if (!isElastic_)
    {
      throw std::runtime_error(
          "DAS receivers require elastic simulation (--is-elastic true)");
    }
    dasGaugeLength_ = opt.das_gauge_length;
    dasDirection_ =
        SourceAndReceiverUtils::ComputeDASVector(opt.das_dip, opt.das_azimuth);
    dasVector_ = dasDirection_;

    // For dipole mode, divide direction by gauge length (GEOS convention)
    if (dasType_ == SourceAndReceiverUtils::DASType::kDipole)
    {
      for (int i = 0; i < 3; ++i) dasVector_[i] /= dasGaugeLength_;
    }

    std::cout << "DAS receiver enabled: type=" << opt.das_type
              << ", dip=" << opt.das_dip << " deg, azimuth=" << opt.das_azimuth
              << " deg, gauge=" << dasGaugeLength_
              << " m, samples=" << dasNumSamples_ << std::endl;
    std::cout << "DAS direction vector: (" << dasDirection_[0] << ", "
              << dasDirection_[1] << ", " << dasDirection_[2] << ")"
              << std::endl;
  }
}

void SEMproxy::run()
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  bool isElastic = isElastic_;

  // Initialize Solver with Partition Info & Compute Local Mass
  m_solver->computeFEInit(*m_mesh, sponge_size_, surface_sponge_, taper_delta_);

  // Synchronize Mass Matrix (Critical for DD)
  if (par_topology_.isDistributed())
  {
    m_syncer->synchronize(m_solver->getMassMatrix(), par_topology_);
    for (int c = 0; c < m_solver->getNumComponents(); ++c)
    {
      m_syncer->synchronize(m_solver->getDampingMatrix(c), par_topology_);
    }
  }

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
        saveSnapshot(indexTimeSample, pnGlobalPrev);

        // Save XY-slice through center Z as plain text for easy visualization
        int nx = nb_nodes_[0];
        int ny = nb_nodes_[1];
        int nz = nb_nodes_[2];
        int zSlice = nz / 2;
        std::ostringstream fname;
        fname << "slice_" << std::setfill('0') << std::setw(5)
              << indexTimeSample << ".dat";
        std::ofstream fslice(fname.str());
        for (int iy = 0; iy < ny; ++iy)
        {
          for (int ix = 0; ix < nx; ++ix)
          {
            int nodeIdx = ix + iy * nx + zSlice * nx * ny;
            fslice << pnGlobalPrev(nodeIdx);
            if (ix < nx - 1) fslice << " ";
          }
          fslice << "\n";
        }
        fslice.close();
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

    // Save receiver trace to plain text file
    {
      std::ofstream fout("receiver_trace.txt");
      fout << "# time pressure_at_receiver\n";
      for (int t = 0; t < num_sample_; ++t)
      {
        fout << t * dt_ << " " << pnAtReceiver(0, t) << "\n";
      }
      fout.close();
      std::cout << "Receiver trace saved to receiver_trace.txt (" << num_sample_
                << " samples)" << std::endl;
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

      // Compute DAS signal at this timestep
      if (dasType_ != SourceAndReceiverUtils::DASType::kNone)
      {
        float dasVal = 0.0f;
        const int totalDASNodes = static_cast<int>(dasNodeIds_.size());
        for (int iNode = 0; iNode < totalDASNodes; ++iNode)
        {
          int nId = dasNodeIds_[iNode];
          if (nId >= 0)
          {
            float w = dasWeights_[iNode];
            dasVal += (uxnGlobalCurr(nId) * dasVector_[0] +
                       uynGlobalCurr(nId) * dasVector_[1] +
                       uznGlobalCurr(nId) * dasVector_[2]) *
                      w;
          }
        }
        dasSignal_(indexTimeSample) = dasVal;
      }

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

    // Save DAS trace to text file
    if (dasType_ != SourceAndReceiverUtils::DASType::kNone)
    {
      std::ofstream fDAS("das_trace.txt");
      fDAS << "# time das_signal\n";
      for (int t = 0; t < num_sample_; ++t)
      {
        fDAS << t * dt_ << " " << dasSignal_(t) << "\n";
      }
      fDAS.close();
      std::cout << "Saved DAS trace to das_trace.txt (" << num_sample_
                << " samples)" << std::endl;
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

  // DAS signal array
  if (dasType_ != SourceAndReceiverUtils::DASType::kNone)
  {
    dasSignal_ = allocateVector<vectorReal>(num_sample_, "dasSignal");
  }
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

  // DAS receiver precomputation
  if (dasType_ != SourceAndReceiverUtils::DASType::kNone)
  {
    const int npe = m_mesh->getNumberOfPointsPerElement();
    const int totalDASNodes = dasNumSamples_ * npe;
    dasNodeIds_.resize(totalDASNodes, -1);
    dasWeights_.resize(totalDASNodes, 0.0f);

    // Compute sample point locations along fiber [-0.5, ..., +0.5]
    std::vector<float> sampleLocations(dasNumSamples_);
    if (dasNumSamples_ == 1)
    {
      sampleLocations[0] = 0.0f;
    }
    else
    {
      for (int i = 0; i < dasNumSamples_; ++i)
      {
        sampleLocations[i] =
            -0.5f + static_cast<float>(i) / (dasNumSamples_ - 1);
      }
    }

    // Compute integration constants
    std::vector<float> integrationConstants(dasNumSamples_);
    if (dasType_ == SourceAndReceiverUtils::DASType::kDipole)
    {
      integrationConstants[0] = -1.0f;
      integrationConstants[1] = 1.0f;
    }
    else
    {
      for (int i = 0; i < dasNumSamples_; ++i)
      {
        integrationConstants[i] = 1.0f / dasNumSamples_;
      }
    }

    // For each sample point along the fiber
    for (int iSample = 0; iSample < dasNumSamples_; ++iSample)
    {
      // Physical coordinates of sample point
      std::array<float, 3> sampleCoord;
      for (int d = 0; d < 3; ++d)
      {
        sampleCoord[d] = rcv_coord_[d] + dasDirection_[d] * dasGaugeLength_ *
                                             sampleLocations[iSample];
      }

      // Find element containing this sample point
      float rel_x = sampleCoord[0] - m_localParams.origin_x;
      int sampleElemIdx =
          static_cast<int>(std::floor((rel_x * ex) / lx)) +
          static_cast<int>(std::floor((sampleCoord[1] * ey) / ly)) * ex +
          static_cast<int>(std::floor((sampleCoord[2] * ez) / lz)) * ey * ex;

      // Clamp to valid range
      if (sampleElemIdx < 0) sampleElemIdx = 0;
      if (sampleElemIdx >= m_mesh->getNumberOfElements()) sampleElemIdx = 0;

      // Get corner coordinates of sample element
      float sampleCornerCoords[8][3];
      int ci = 0;
      for (int kk : nodes_corner)
      {
        for (int jj : nodes_corner)
        {
          for (int ii : nodes_corner)
          {
            int nIdx = m_mesh->globalNodeIndex(sampleElemIdx, ii, jj, kk);
            sampleCornerCoords[ci][0] = m_mesh->nodeCoord(nIdx, 0);
            sampleCornerCoords[ci][1] = m_mesh->nodeCoord(nIdx, 1);
            sampleCornerCoords[ci][2] = m_mesh->nodeCoord(nIdx, 2);
            ci++;
          }
        }
      }

      // Store global node IDs for this sample's element
      int baseIdx = iSample * npe;
      for (int kk = 0; kk < order + 1; ++kk)
      {
        for (int jj = 0; jj < order + 1; ++jj)
        {
          for (int ii = 0; ii < order + 1; ++ii)
          {
            int localNode =
                ii + jj * (order + 1) + kk * (order + 1) * (order + 1);
            dasNodeIds_[baseIdx + localNode] =
                m_mesh->globalNodeIndex(sampleElemIdx, ii, jj, kk);
          }
        }
      }

      // Compute weights for this sample point
      switch (order)
      {
        case 1:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<1>(
              sampleCornerCoords, sampleCoord, dasDirection_,
              integrationConstants[iSample], dasType_, &dasWeights_[baseIdx]);
          break;
        case 2:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<2>(
              sampleCornerCoords, sampleCoord, dasDirection_,
              integrationConstants[iSample], dasType_, &dasWeights_[baseIdx]);
          break;
        case 3:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<3>(
              sampleCornerCoords, sampleCoord, dasDirection_,
              integrationConstants[iSample], dasType_, &dasWeights_[baseIdx]);
          break;
        default:
          throw std::runtime_error("Unsupported order for DAS: " +
                                   std::to_string(order));
      }
    }

    std::cout << "DAS precomputation complete: " << dasNumSamples_
              << " samples, " << totalDASNodes << " node entries" << std::endl;
  }
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
