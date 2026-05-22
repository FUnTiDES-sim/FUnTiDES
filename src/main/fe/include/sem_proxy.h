/**
 * @file sem_proxy.h
 * @brief Defines the main interface and orchestrator for the Spectral Element Method (SEM) proxy simulation.
 */

#ifndef FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
#define FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_

#include <data_type.h>
#include <utils.h>

#include <array>
#include <chrono>
#include <future>
#include <memory>
#include <string>
#include <vector>

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
 *
 * Orchestrates the entire simulation lifecycle, including MPI setup,
 * mesh parameters, memory allocation (Device/Host), solver initialization,
 * time-stepping, and asynchronous I/O operations.
 */
class SEMproxy {
 public:
  /**
   * @brief Constructs the SEMproxy simulation environment.
   * @param cfg The parsed configuration options for the simulation.
   */
  explicit SEMproxy(const SemProxyOptions& cfg);

  /**
   * @brief Destroys the SEMproxy instance, safely releasing I/O controllers.
   */
  ~SEMproxy() { io_ctrl_.reset(); }

  /**
   * @brief Initializes the finite element arrays, host mirrors, and compute sources.
   */
  void InitFiniteElem() {
    InitArrays();
    InitSource();
  };

  /**
   * @brief Executes the main time-stepping compute loop of the simulation.
   */
  void Run();

  /**
   * @brief Saves a 2D slice of the computational domain to a file.
   * @param host_slice The host-side data array containing the slice data.
   * @param size_x The number of elements along the X axis.
   * @param size_y The number of elements along the Y axis.
   * @param filepath The destination file path.
   */
  void SaveSlice(const vectorReal& host_slice, int size_x, int size_y, const std::string& filepath) const;

  /**
   * @brief Dispatches a deep copy from Device to Host and asynchronously saves the snapshot.
   * @param time_sample The current iteration step.
   * @param d_data The Device view containing the data to save.
   * @param h_data The pre-allocated Host mirror to copy data into.
   */
  void SaveSnapshot(int time_sample, const vectorReal& d_data, vectorReal::host_mirror_type& h_data) const;

  /**
   * @brief Computes a stable time step (dt) automatically based on the CFL condition.
   * @param cfl_factor The Courant-Friedrichs-Lewy stability factor.
   * @return The computed time step in seconds.
   */
  float FindCflDt(float cfl_factor);

 private:
  model::CartesianParams<float, int> local_params_;  ///< Local Cartesian parameters for the subdomain.
  utils::DistributedContext dist_ctx_;               ///< MPI distributed context (rank, size).
  utils::ParallelTopology par_topology_;             ///< Parallel topology layout.

  int num_elements_[3] = {0};   ///< Number of elements along [X, Y, Z] for this MPI rank.
  int num_nodes_[3] = {0};      ///< Number of nodes along [X, Y, Z] for this MPI rank.
  float domain_size_[3] = {0};  ///< Physical dimensions of the local domain.

  // --- I/O ---
  bool is_snapshots_ = false;                 ///< Flag indicating if 3D snapshots are enabled.
  int snap_time_interval_ = 0;                ///< Number of iterations between snapshots.
  std::string snap_folder_;                   ///< Directory path for saving snapshots.
  std::shared_ptr<SemIOController> io_ctrl_;  ///< Controller for ADIOS2 I/O operations.

  // --- ASYNC I/O ---
  std::vector<std::future<void>> snapshot_futures_;  ///< Tracks background async I/O threads.
  /**
   * @brief Blocks the main thread until all pending async I/O tasks complete.
   */
  void WaitSnapshots();

  // --- Physics & Meshing Flags ---
  bool is_elastic_ = false;          ///< True if simulating elastic wave propagation.
  bool is_acousto_elastic_ = false;  ///< True if simulating coupled acousto-elastic wave propagation.
  bool free_surface_ = false;        ///< True if the top boundary acts as a free surface.
  bool is_dg_ = false;               ///< True if using Discontinuous Galerkin method.
  bool is_dg_sem_ = false;           ///< True if using Discontinuous Galerkin - Spectral Element method coupling.

  std::array<float, 3> sponge_size_ = {0, 0, 0};  ///< Thickness of absorbing boundaries (sponge layers).
  bool surface_sponge_ = false;                   ///< True if the top surface has an absorbing boundary.
  float taper_delta_ = 0.015f;                    ///< Tapering coefficient for the sponge boundaries.

  float dt_ = 0.0f;        ///< Time step size (seconds).
  float time_max_ = 0.0f;  ///< Maximum simulation time (seconds).
  int num_samples_ = 0;    ///< Total number of time steps to compute.

  const int num_rhs_ = 1;   ///< Number of Right-Hand Side source terms.
  int source_element_ = 0;  ///< Index of the element containing the source.
  float t_peak_ = 0.0f;     ///< Time peak of the Ricker wavelet source.
  float f0_ = 0.0f;         ///< Dominant frequency of the source in Hz.
  int ricker_order_ = 0;    ///< Order of the Ricker wavelet.

  std::array<float, 3> src_coord_ = {0};  ///< Global coordinates of the source (X, Y, Z).
  std::array<float, 3> rcv_coord_ = {0};  ///< Global coordinates of the receiver (X, Y, Z).

  std::shared_ptr<model::ModelApi<float, int>> mesh_;         ///< Pointer to the finite element mesh API.
  std::unique_ptr<solver::fe::Solver> solver_;                ///< Main numerical solver instance.
  std::unique_ptr<solver::fe::BoundarySynchronizer> syncer_;  ///< Handles MPI boundary ghost-node synchronization.
  SolverUtils utils_;                                         ///< General solver math utilities.

  // --- Acoustic / Shared Arrays (Device) ---
  arrayReal rhs_term_;         ///< Source term array over time (Device).
  arrayReal rhs_term_dg_;      ///< DG source term array over time (Device).
  arrayReal rhs_term_sem_;     ///< SEM source term array over time (Device).
  vectorReal pn_global_prev_;  ///< Pressure field at time t-1 (Device).
  vectorReal pn_global_curr_;  ///< Pressure field at time t (Device).
  arrayReal pn_dg_prev_;       ///< DG Pressure field at time t-1 (Device).
  arrayReal pn_dg_curr_;       ///< DG Pressure field at time t (Device).
  vectorReal pn_sem_prev_;     ///< SEM Pressure field at time t-1 (Device).
  vectorReal pn_sem_curr_;     ///< SEM Pressure field at time t (Device).
  vectorInt rhs_element_;      ///< Element indices containing sources (Device).
  vectorInt rhs_element_rcv_;  ///< Element indices containing receivers (Device).
  arrayReal rhs_weights_;      ///< Interpolation weights for sources (Device).
  arrayReal rhs_weights_rcv_;  ///< Interpolation weights for receivers (Device).
  arrayReal pn_at_receiver_;   ///< Recorded pressure traces at receivers (Device).

  // --- Elastic Arrays (Device) ---
  arrayReal rhs_term_x_;        ///< X-component of the source term (Device).
  arrayReal rhs_term_y_;        ///< Y-component of the source term (Device).
  arrayReal rhs_term_z_;        ///< Z-component of the source term (Device).
  vectorReal uxn_global_prev_;  ///< X-displacement at time t-1 (Device).
  vectorReal uyn_global_prev_;  ///< Y-displacement at time t-1 (Device).
  vectorReal uzn_global_prev_;  ///< Z-displacement at time t-1 (Device).
  vectorReal uxn_global_curr_;  ///< X-displacement at time t (Device).
  vectorReal uyn_global_curr_;  ///< Y-displacement at time t (Device).
  vectorReal uzn_global_curr_;  ///< Z-displacement at time t (Device).
  arrayReal uxn_at_receiver_;   ///< Recorded X-displacement traces at receivers (Device).
  arrayReal uyn_at_receiver_;   ///< Recorded Y-displacement traces at receivers (Device).
  arrayReal uzn_at_receiver_;   ///< Recorded Z-displacement traces at receivers (Device).

  // --- DAS Receiver Data ---
  SourceAndReceiverUtils::DASType das_type_ = SourceAndReceiverUtils::DASType::kNone;  ///< Type of DAS receiver.
  int das_num_samples_ = 5;                         ///< Number of integration samples along the fiber.
  float das_gauge_length_ = 1.0f;                   ///< Gauge length for the DAS fiber in meters.
  std::array<float, 3> das_direction_ = {1, 0, 0};  ///< Fiber direction unit vector.
  std::array<float, 3> das_vector_ = {1, 0, 0};     ///< Scaled fiber direction vector.
  std::vector<int> das_node_ids_;                   ///< Node IDs involved in DAS integration.
  std::vector<float> das_weights_;                  ///< Weights for DAS integration points.
  vectorReal das_signal_;                           ///< Output DAS signal trace over time (Device).

  // --- HOST MIRRORS (Used to avoid UVM overhead when CPU needs data) ---
  vectorInt::host_mirror_type h_rhs_element_;      ///< CPU mirror for source elements.
  vectorInt::host_mirror_type h_rhs_element_rcv_;  ///< CPU mirror for receiver elements.
  arrayReal::host_mirror_type h_rhs_weights_;      ///< CPU mirror for source interpolation weights.
  arrayReal::host_mirror_type h_rhs_weights_rcv_;  ///< CPU mirror for receiver interpolation weights.
  arrayReal::host_mirror_type h_rhs_term_;         ///< CPU mirror for acoustic source term.
  arrayReal::host_mirror_type h_rhs_term_dg_;      ///< CPU mirror for DG-SEM DG source term.
  arrayReal::host_mirror_type h_rhs_term_sem_;     ///< CPU mirror for DG-SEM SEM source term.
  arrayReal::host_mirror_type h_rhs_term_x_;       ///< CPU mirror for elastic X source term.
  arrayReal::host_mirror_type h_rhs_term_y_;       ///< CPU mirror for elastic Y source term.
  arrayReal::host_mirror_type h_rhs_term_z_;       ///< CPU mirror for elastic Z source term.

  arrayReal::host_mirror_type h_pn_at_receiver_;   ///< CPU mirror for acoustic receiver traces.
  arrayReal::host_mirror_type h_uxn_at_receiver_;  ///< CPU mirror for elastic X receiver traces.
  arrayReal::host_mirror_type h_uyn_at_receiver_;  ///< CPU mirror for elastic Y receiver traces.
  arrayReal::host_mirror_type h_uzn_at_receiver_;  ///< CPU mirror for elastic Z receiver traces.
  vectorReal::host_mirror_type h_das_signal_;      ///< CPU mirror for DAS signal trace.

  vectorReal::host_mirror_type h_pn_global_curr_;   ///< CPU mirror for current pressure field.
  vectorReal::host_mirror_type h_pn_global_prev_;   ///< CPU mirror for previous pressure field.
  vectorReal::host_mirror_type h_pn_sem_curr_;      ///< CPU mirror for SEM current pressure field.
  vectorReal::host_mirror_type h_pn_sem_prev_;      ///< CPU mirror for SEM previous pressure field.
  arrayReal::host_mirror_type h_pn_dg_curr_;        ///< CPU mirror for DG current pressure field.
  arrayReal::host_mirror_type h_pn_dg_prev_;        ///< CPU mirror for DG previous pressure field.
  vectorReal::host_mirror_type h_uxn_global_curr_;  ///< CPU mirror for current X-displacement.
  vectorReal::host_mirror_type h_uyn_global_curr_;  ///< CPU mirror for current Y-displacement.
  vectorReal::host_mirror_type h_uzn_global_curr_;  ///< CPU mirror for current Z-displacement.
  vectorReal::host_mirror_type h_uxn_global_prev_;  ///< CPU mirror for previous X-displacement.
  vectorReal::host_mirror_type h_uyn_global_prev_;  ///< CPU mirror for previous Y-displacement.
  vectorReal::host_mirror_type h_uzn_global_prev_;  ///< CPU mirror for previous Z-displacement.

  // --- Performance Tracking ---
  double time_init_ = 0.0;     ///< Total initialization time in seconds.
  double time_compute_ = 0.0;  ///< Total solver execution time in seconds.
  double time_io_ = 0.0;       ///< Total I/O operations time in seconds.

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

  // --- Options Translators ---
  int GetPhysic(std::string physic_arg);
  utils::enums::implemType GetImplem(std::string implem_arg);
  utils::enums::methodType GetMethod(std::string method_arg);
  utils::enums::meshType GetMesh(std::string mesh_arg);
  model::AnisotropyType GetAnisotropy(std::string anisotropy_arg);
};

#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
