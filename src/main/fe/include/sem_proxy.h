// This file defines the main interface for the SEM proxy simulation.

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

class SEMproxy {
 public:
  explicit SEMproxy(const SemProxyOptions& cfg);

  ~SEMproxy() { io_ctrl_.reset(); }

  void InitFiniteElem() {
    InitArrays();
    InitSource();
  };

  void Run();

  void SaveSlice(const vectorReal& host_slice, int size_x, int size_y, const std::string& filepath) const;
  void SaveSnapshot(int time_sample, const vectorReal& d_data, vectorReal::host_mirror_type& h_data) const;
  float FindCflDt(float cfl_factor);

 private:
  model::CartesianParams<float, int> local_params_;
  utils::DistributedContext dist_ctx_;
  utils::ParallelTopology par_topology_;

  int num_elements_[3] = {0};
  int num_nodes_[3] = {0};
  float domain_size_[3] = {0};

  // --- I/O ---
  bool is_snapshots_ = false;
  int snap_time_interval_ = 0;
  std::string snap_folder_;
  std::shared_ptr<SemIOController> io_ctrl_;
  // --- ASYNC I/O ---
  std::vector<std::future<void>> snapshot_futures_;
  void WaitSnapshots();

  bool is_elastic_ = false;
  bool is_acousto_elastic_ = false;
  bool free_surface_ = false;
  bool is_dg_ = false;

  std::array<float, 3> sponge_size_ = {0, 0, 0};
  bool surface_sponge_ = false;
  float taper_delta_ = 0.015f;

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

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<solver::fe::Solver> solver_;
  std::unique_ptr<solver::fe::BoundarySynchronizer> syncer_;
  SolverUtils utils_;

  // --- Acoustic / Shared Arrays (Device) ---
  arrayReal rhs_term_;
  vectorReal pn_global_prev_;
  vectorReal pn_global_curr_;
  arrayReal pn_dg_prev_;
  arrayReal pn_dg_curr_;
  vectorInt rhs_element_;
  vectorInt rhs_element_rcv_;
  arrayReal rhs_weights_;
  arrayReal rhs_weights_rcv_;
  arrayReal pn_at_receiver_;

  // --- Elastic Arrays (Device) ---
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

  // --- HOST MIRRORS (Used to avoid UVM overhead when CPU needs data) ---
  vectorInt::host_mirror_type h_rhs_element_;
  vectorInt::host_mirror_type h_rhs_element_rcv_;
  arrayReal::host_mirror_type h_rhs_weights_;
  arrayReal::host_mirror_type h_rhs_weights_rcv_;
  arrayReal::host_mirror_type h_rhs_term_;
  arrayReal::host_mirror_type h_rhs_term_x_;
  arrayReal::host_mirror_type h_rhs_term_y_;
  arrayReal::host_mirror_type h_rhs_term_z_;

  arrayReal::host_mirror_type h_pn_at_receiver_;
  arrayReal::host_mirror_type h_uxn_at_receiver_;
  arrayReal::host_mirror_type h_uyn_at_receiver_;
  arrayReal::host_mirror_type h_uzn_at_receiver_;
  vectorReal::host_mirror_type h_das_signal_;

  vectorReal::host_mirror_type h_pn_global_curr_;
  vectorReal::host_mirror_type h_pn_global_prev_;
  vectorReal::host_mirror_type h_uxn_global_curr_;
  vectorReal::host_mirror_type h_uyn_global_curr_;
  vectorReal::host_mirror_type h_uzn_global_curr_;
  vectorReal::host_mirror_type h_uxn_global_prev_;
  vectorReal::host_mirror_type h_uyn_global_prev_;
  vectorReal::host_mirror_type h_uzn_global_prev_;

  // --- Performance Tracking ---
  double time_init_ = 0.0;
  double time_compute_ = 0.0;
  double time_io_ = 0.0;

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

  int GetPhysic(std::string physic_arg);
  utils::enums::implemType GetImplem(std::string implem_arg);
  utils::enums::methodType GetMethod(std::string method_arg);
  utils::enums::meshType GetMesh(std::string mesh_arg);
  model::AnisotropyType GetAnisotropy(std::string anisotropy_arg);
};

#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
