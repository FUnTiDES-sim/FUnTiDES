// This file defines the main interface for the SEM proxy simulation.

#ifndef FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
#define FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_

#include <data_type.h>
#include <utils.h>

#include <array>
#include <chrono>
#include <memory>
#include <string>

#include "boundary_synchronizer.h"
#include "cartesian_params.h"
#include "distributed_ctx.h"
#include "model_struct.h"
#include "model_unstruct.h"
#include "sem_enums.h"
#include "sem_io_controller.h"
#include "sem_proxy_options.h"
#include "sem_solver.h"
#include "solver_factory.h"
#include "source_and_receiver_utils.h"

/**
 * @class SEMproxy
 * @brief Main driver class for the Spectral Element Method (SEM) simulation.
 */
class SEMproxy {
 public:
  /**
   * @brief Constructs the SEMproxy simulation environment.
   * @param cfg The parsed configuration options for the simulation.
   */
  explicit SEMproxy(const SemProxyOptions& cfg);

  /**
   * @brief Destructs the SEMproxy instance and closes internal engines.
   */
  ~SEMproxy() { io_ctrl_.reset(); }

  /**
   * @brief Initializes the finite element arrays and sources.
   */
  void InitFiniteElem() {
    InitArrays();
    InitSource();
  };

  /**
   * @brief Executes the main time-stepping simulation loop.
   */
  void Run();

  /**
   * @brief Saves a 2D data slice in Gnuplot matrix format.
   * @param host_slice The data view to save.
   * @param size_x X-dimension size.
   * @param size_y Y-dimension size.
   * @param filepath The output file path.
   */
  void SaveSlice(const vectorReal& host_slice, int size_x, int size_y, const std::string& filepath) const;

  /**
   * @brief Saves a full simulation snapshot for a given time sample.
   * @param time_sample The current timestep index.
   * @param data The data array representing the field state.
   */
  void SaveSnapshot(int time_sample, vectorReal data) const;
  void SaveSnapshot(int timestep);

  /**
   * @brief Computes optimal time step using the CFL stability condition.
   * @param cfl_factor Stability factor (e.g., 0.5-0.7 for 2nd-order schemes).
   * @return The maximum stable time step delta (dt).
   */
  float FindCflDt(float cfl_factor);

 private:
  // --- Domain Decomposition & Topology ---
  model::CartesianParams<float, int> local_params_;
  utils::DistributedContext dist_ctx_;
  utils::ParallelTopology par_topology_;

  int num_elements_[3] = {0};
  int num_nodes_[3] = {0};
  float domain_size_[3] = {0};

  // --- IO & Snapshots ---
  bool is_snapshots_ = false;
  int snap_time_interval_ = 0;
  std::string snap_folder_;
  std::shared_ptr<SemIOController> io_ctrl_;

  // --- Physics Configurations ---
  bool is_elastic_ = false;
  bool is_acousto_elastic_ = false;
  bool free_surface_ = false;

  // --- Sponge Boundary Parameters ---
  std::array<float, 3> sponge_size_ = {0, 0, 0};
  bool surface_sponge_ = false;
  float taper_delta_ = 0.015f;

  // --- Time & Source Parameters ---
  float dt_ = 0.0f;
  float time_max_ = 0.0f;
  int num_samples_ = 0;

  const int num_rhs_ = 1;
  int source_element_ = 0;
  float t_peak_ = 0.0f;
  float f0_ = 0.0f;
  int ricker_order_ = 0;

  std::array<float, 3> src_coord_ = {0};
  std::array<float, 3> rcv_coord_ = {0};

  // --- Core Solver & Mesh ---
  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<solver::fe::Solver> solver_;
  std::unique_ptr<solver::fe::BoundarySynchronizer> syncer_;
  SolverUtils utils_;

  // --- Acoustic / Shared Arrays ---
  arrayReal rhs_term_;
  vectorReal pn_global_prev_;
  vectorReal pn_global_curr_;
  vectorInt rhs_element_;
  vectorInt rhs_element_rcv_;
  arrayReal rhs_weights_;
  arrayReal rhs_weights_rcv_;
  arrayReal pn_at_receiver_;

  // --- Elastic Arrays ---
  arrayReal rhs_term_x_;
  arrayReal rhs_term_y_;
  arrayReal rhs_term_z_;
  vectorReal uxn_global_prev_;
  vectorReal uyn_global_prev_;
  vectorReal uzn_global_prev_;
  vectorReal uxn_global_curr_;
  vectorReal uyn_global_curr_;
  vectorReal uzn_global_curr_;
  arrayReal uxn_at_receiver_;
  arrayReal uyn_at_receiver_;
  arrayReal uzn_at_receiver_;

  // --- DAS Receiver Data ---
  SourceAndReceiverUtils::DASType das_type_ = SourceAndReceiverUtils::DASType::kNone;
  int das_num_samples_ = 5;
  float das_gauge_length_ = 1.0f;
  std::array<float, 3> das_direction_ = {1, 0, 0};
  std::array<float, 3> das_vector_ = {1, 0, 0};
  std::vector<int> das_node_ids_;
  std::vector<float> das_weights_;
  vectorReal das_signal_;

  // --- Performance Tracking ---
  double time_init_ = 0.0;
  double time_compute_ = 0.0;
  double time_io_ = 0.0;

  // --- Initialization Helpers ---
  void InitSource();
  void InitArrays();
  void InitMpi(int* mpi_init);
  void InitSimParams(const SemProxyOptions& opt);
  void InitMeshParams(const SemProxyOptions& opt);
  void InitTopology();
  void InitSync();
  void InitTimeParams(const SemProxyOptions& opt);

  void SetupSolver(const SemProxyOptions& opt);
  void SetupAttenuation(const SemProxyOptions& opt);
  void SetupIo(const SemProxyOptions& opt);
  void SetupDas(const SemProxyOptions& opt);

  void DisplayInitMsg(const SemProxyOptions& opt);
  void DisplayPerfMsg() const;

  // --- Parsing Helpers ---
  int GetPhysic(std::string physic_arg);
  utils::enums::implemType GetImplem(std::string implem_arg);
  utils::enums::methodType GetMethod(std::string method_arg);
  utils::enums::meshType GetMesh(std::string mesh_arg);
  model::AnisotropyType GetAnisotropy(std::string anisotropy_arg);
};

#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
