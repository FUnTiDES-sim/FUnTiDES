/**
 * @file sem_proxy.h
 * @brief Driver of a SEM/DG wave-propagation run: setup, time loop, I/O and timing.
 */

#ifndef FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
#define FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_

#include <data_type.h>
#include <source_time_function.h>

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
 * @brief Owns the mesh, the solver, the source/receiver data and the device
 *        fields of one simulation, and runs the time loop.
 *
 * Built from a parsed option set. After construction, InitFiniteElem() must be
 * called before Run(). Snapshots are written asynchronously; pending writes are
 * awaited internally.
 */
class SEMproxy {
 public:
  /**
   * @brief Builds the simulation from the parsed options.
   * @param[in] cfg Parsed configuration options.
   */
  explicit SEMproxy(const SemProxyOptions& cfg);

  /**
   * @brief Releases the I/O controller.
   */
  ~SEMproxy() { io_ctrl_.reset(); }

  /**
   * @brief Allocates the solution arrays and host mirrors and computes the source terms.
   */
  void InitFiniteElem() {
    InitArrays();
    InitSource();
  };

  /**
   * @brief Runs the time loop for num_samples_ steps.
   */
  void Run();

  /**
   * @brief Writes a 2D slice of the domain to a file.
   * @param[in] host_slice Slice values on the host.
   * @param[in] size_x Extent of the slice along X.
   * @param[in] size_y Extent of the slice along Y.
   * @param[in] filepath Destination file.
   * @todo VERIFY: are size_x and size_y counted in nodes or elements, and what is the storage order and file format of
   * host_slice?
   */
  void SaveSlice(const vectorReal& host_slice, int size_x, int size_y, const std::string& filepath) const;

  /**
   * @brief Copies a device field to the host and writes it as a snapshot in a background task.
   * @param[in] time_sample Index of the current time step.
   * @param[in] d_data Device field to save.
   * @param[in,out] h_data Preallocated host mirror used as the copy destination.
   */
  void SaveSnapshot(int time_sample, const vectorReal& d_data, vectorReal::host_mirror_type& h_data) const;

  /**
   * @brief Computes a time step from a CFL stability factor.
   * @param[in] cfl_factor Courant-Friedrichs-Lewy factor.
   * @return Time step in seconds.
   * @todo VERIFY: which spacing and speed (min spacing, max speed) enter the estimate?
   */
  float FindCflDt(float cfl_factor);

 private:
  model::CartesianParams<float, int> local_params_;  ///< Cartesian parameters of the local subdomain.
  utils::DistributedContext dist_ctx_;               ///< MPI rank and size.
  utils::ParallelTopology par_topology_;             ///< Layout of the ranks.

  int num_elements_[3] = {0};   ///< Number of elements along x, y, z on this rank.
  int num_nodes_[3] = {0};      ///< Number of nodes along x, y, z on this rank.
  float domain_size_[3] = {0};  ///< Extent of the local domain along x, y, z. @todo VERIFY: unit (meters?).

  // Snapshot I/O.
  bool is_snapshots_ = false;                 ///< True if 3D snapshots are written.
  int snap_time_interval_ = 0;                ///< Number of time steps between two snapshots.
  std::string snap_folder_;                   ///< Output directory of the snapshots.
  std::shared_ptr<SemIOController> io_ctrl_;  ///< ADIOS2 controller used for snapshots and receivers.

  // Asynchronous I/O.
  std::vector<std::future<void>> snapshot_futures_;  ///< Pending background snapshot writes.

  /**
   * @brief Blocks until all pending snapshot writes are complete.
   */
  void WaitSnapshots();

  // Physics and method selection.
  bool is_elastic_ = false;              ///< True for elastic propagation.
  bool is_acousto_elastic_ = false;      ///< True for coupled acousto-elastic propagation.
  bool free_surface_ = false;            ///< True if the top boundary is a free surface.
  bool is_dg_ = false;                   ///< True if the DG method is used.
  bool is_dg_sem_ = false;               ///< True if DG is coupled with SEM.
  bool is_dg_padaptive_ = false;         ///< True if the p-adaptive DG method is used.
  float dg_sem_iface_z_ = 1000.f;        ///< z coordinate of the DG-SEM interface.
  float dg_padaptive_iface_z_ = 1000.f;  ///< z coordinate of the pMin/pMax interface of the p-adaptive DG method.
  int order_min_ = 0;                    ///< Lower polynomial order of the p-adaptive DG method.

  std::array<float, 3> sponge_size_ = {0, 0, 0};  ///< Thickness of the sponge layers along x, y, z.
  bool surface_sponge_ = false;  ///< See Solver::computeFEInit. @todo VERIFY: is the top surface absorbing when true?
  float taper_delta_ = 0.015f;   ///< Taper coefficient of the sponge damping.

  float dt_ = 0.0f;        ///< Time step in seconds.
  float time_max_ = 0.0f;  ///< Simulated duration in seconds.
  int num_samples_ = 0;    ///< Number of time steps.

  const int num_rhs_ = 1;   ///< Number of sources.
  int source_element_ = 0;  ///< Index of the element containing the source.
  float t_peak_ =
      0.0f;  ///< Time of the peak of the source wavelet. @todo VERIFY: unit and reference (seconds from t = 0?).
  float f0_ = 0.0f;       ///< Dominant frequency of the source in Hz.
  int ricker_order_ = 0;  ///< Derivative order of the Ricker wavelet.

  std::array<float, 3> src_coord_ = {0};  ///< Global coordinates (x, y, z) of the source.
  std::array<float, 3> rcv_coord_ = {0};  ///< Global coordinates (x, y, z) of the receiver.

  std::shared_ptr<model::ModelApi<float, int>> mesh_;         ///< Mesh and model of the local subdomain.
  std::unique_ptr<solver::fe::Solver> solver_;                ///< Solver advancing the fields.
  std::unique_ptr<solver::fe::BoundarySynchronizer> syncer_;  ///< Exchanges boundary nodes between MPI ranks.
  SourceTimeFunction source_time_function_;                   ///< Source wavelet.

  // Acoustic and shared arrays (device).
  arrayReal rhs_term_;              ///< Source term over time.
  arrayReal rhs_term_dg_;           ///< Source term over time, DG part.
  arrayReal rhs_term_sem_;          ///< Source term over time, SEM part.
  arrayReal rhs_term_pmin_;         ///< Source term over time, pMin domain.
  arrayReal rhs_term_pmax_;         ///< Source term over time, pMax domain.
  vectorReal pn_global_prev_;       ///< Pressure at time step n-1.
  vectorReal pn_global_curr_;       ///< Pressure at time step n.
  arrayReal pn_dg_prev_;            ///< DG pressure at time step n-1.
  arrayReal pn_dg_curr_;            ///< DG pressure at time step n.
  vectorReal pn_sem_prev_;          ///< SEM pressure at time step n-1.
  vectorReal pn_sem_curr_;          ///< SEM pressure at time step n.
  arrayReal pn_pmin_dg_prev_;       ///< pMin DG pressure at time step n-1.
  arrayReal pn_pmin_dg_curr_;       ///< pMin DG pressure at time step n.
  arrayReal pn_pmax_dg_prev_;       ///< pMax DG pressure at time step n-1.
  arrayReal pn_pmax_dg_curr_;       ///< pMax DG pressure at time step n.
  vectorInt rhs_element_;           ///< Elements containing the sources.
  vectorInt rhs_element_rcv_;       ///< Elements containing the receivers.
  arrayReal rhs_weights_;           ///< Interpolation weights of the sources.
  arrayReal rhs_weights_rcv_;       ///< Interpolation weights of the receivers.
  arrayReal rhs_pmin_weights_;      ///< Source interpolation weights, pMin domain.
  arrayReal rhs_pmax_weights_;      ///< Source interpolation weights, pMax domain.
  arrayReal rhs_pmin_weights_rcv_;  ///< Receiver interpolation weights, pMin domain.
  arrayReal rhs_pmax_weights_rcv_;  ///< Receiver interpolation weights, pMax domain.
  arrayReal pn_at_receiver_;        ///< Pressure traces recorded at the receivers.

  // Elastic arrays (device).
  arrayReal rhs_term_x_;        ///< X component of the source term.
  arrayReal rhs_term_y_;        ///< Y component of the source term.
  arrayReal rhs_term_z_;        ///< Z component of the source term.
  vectorReal uxn_global_prev_;  ///< X displacement at time step n-1.
  vectorReal uyn_global_prev_;  ///< Y displacement at time step n-1.
  vectorReal uzn_global_prev_;  ///< Z displacement at time step n-1.
  vectorReal uxn_global_curr_;  ///< X displacement at time step n.
  vectorReal uyn_global_curr_;  ///< Y displacement at time step n.
  vectorReal uzn_global_curr_;  ///< Z displacement at time step n.
  arrayReal uxn_at_receiver_;   ///< X displacement traces recorded at the receivers.
  arrayReal uyn_at_receiver_;   ///< Y displacement traces recorded at the receivers.
  arrayReal uzn_at_receiver_;   ///< Z displacement traces recorded at the receivers.

  // DAS receiver.
  SourceAndReceiverUtils::DASType das_type_ = SourceAndReceiverUtils::DASType::kNone;  ///< DAS receiver type.
  int das_num_samples_ = 5;                         ///< Number of integration samples along the fiber.
  float das_gauge_length_ = 1.0f;                   ///< Gauge length of the fiber in meters.
  std::array<float, 3> das_direction_ = {1, 0, 0};  ///< Unit vector along the fiber.
  std::array<float, 3> das_vector_ = {
      1, 0, 0};                     ///< Fiber direction scaled by a length. @todo VERIFY: which length (gauge length?).
  std::vector<int> das_node_ids_;   ///< Global node indices used by the DAS integration.
  std::vector<float> das_weights_;  ///< Weights of the DAS integration points.
  vectorReal das_signal_;           ///< DAS signal over time (device).

  // Host mirrors of the device arrays above, used when the CPU reads or writes data.
  vectorInt::host_mirror_type h_rhs_element_;           ///< Mirror of rhs_element_.
  vectorInt::host_mirror_type h_rhs_element_rcv_;       ///< Mirror of rhs_element_rcv_.
  arrayReal::host_mirror_type h_rhs_weights_;           ///< Mirror of rhs_weights_.
  arrayReal::host_mirror_type h_rhs_weights_rcv_;       ///< Mirror of rhs_weights_rcv_.
  arrayReal::host_mirror_type h_rhs_term_;              ///< Mirror of rhs_term_.
  arrayReal::host_mirror_type h_rhs_term_dg_;           ///< Mirror of rhs_term_dg_.
  arrayReal::host_mirror_type h_rhs_term_sem_;          ///< Mirror of rhs_term_sem_.
  arrayReal::host_mirror_type h_rhs_term_pmin_;         ///< Mirror of rhs_term_pmin_.
  arrayReal::host_mirror_type h_rhs_term_pmax_;         ///< Mirror of rhs_term_pmax_.
  arrayReal::host_mirror_type h_rhs_pmin_weights_;      ///< Mirror of rhs_pmin_weights_.
  arrayReal::host_mirror_type h_rhs_pmax_weights_;      ///< Mirror of rhs_pmax_weights_.
  arrayReal::host_mirror_type h_rhs_pmin_weights_rcv_;  ///< Mirror of rhs_pmin_weights_rcv_.
  arrayReal::host_mirror_type h_rhs_pmax_weights_rcv_;  ///< Mirror of rhs_pmax_weights_rcv_.
  arrayReal::host_mirror_type h_rhs_term_x_;            ///< Mirror of rhs_term_x_.
  arrayReal::host_mirror_type h_rhs_term_y_;            ///< Mirror of rhs_term_y_.
  arrayReal::host_mirror_type h_rhs_term_z_;            ///< Mirror of rhs_term_z_.

  arrayReal::host_mirror_type h_pn_at_receiver_;   ///< Mirror of pn_at_receiver_.
  arrayReal::host_mirror_type h_uxn_at_receiver_;  ///< Mirror of uxn_at_receiver_.
  arrayReal::host_mirror_type h_uyn_at_receiver_;  ///< Mirror of uyn_at_receiver_.
  arrayReal::host_mirror_type h_uzn_at_receiver_;  ///< Mirror of uzn_at_receiver_.
  vectorReal::host_mirror_type h_das_signal_;      ///< Mirror of das_signal_.

  vectorReal::host_mirror_type h_pn_global_curr_;   ///< Mirror of pn_global_curr_.
  vectorReal::host_mirror_type h_pn_global_prev_;   ///< Mirror of pn_global_prev_.
  vectorReal::host_mirror_type h_pn_sem_curr_;      ///< Mirror of pn_sem_curr_.
  vectorReal::host_mirror_type h_pn_sem_prev_;      ///< Mirror of pn_sem_prev_.
  arrayReal::host_mirror_type h_pn_dg_curr_;        ///< Mirror of pn_dg_curr_.
  arrayReal::host_mirror_type h_pn_dg_prev_;        ///< Mirror of pn_dg_prev_.
  arrayReal::host_mirror_type h_pn_pmin_dg_curr_;   ///< Mirror of pn_pmin_dg_curr_.
  arrayReal::host_mirror_type h_pn_pmin_dg_prev_;   ///< Mirror of pn_pmin_dg_prev_.
  arrayReal::host_mirror_type h_pn_pmax_dg_curr_;   ///< Mirror of pn_pmax_dg_curr_.
  arrayReal::host_mirror_type h_pn_pmax_dg_prev_;   ///< Mirror of pn_pmax_dg_prev_.
  vectorReal::host_mirror_type h_uxn_global_curr_;  ///< Mirror of uxn_global_curr_.
  vectorReal::host_mirror_type h_uyn_global_curr_;  ///< Mirror of uyn_global_curr_.
  vectorReal::host_mirror_type h_uzn_global_curr_;  ///< Mirror of uzn_global_curr_.
  vectorReal::host_mirror_type h_uxn_global_prev_;  ///< Mirror of uxn_global_prev_.
  vectorReal::host_mirror_type h_uyn_global_prev_;  ///< Mirror of uyn_global_prev_.
  vectorReal::host_mirror_type h_uzn_global_prev_;  ///< Mirror of uzn_global_prev_.

  // Timing.
  double time_init_ = 0.0;     ///< Initialization time in seconds.
  double time_compute_ = 0.0;  ///< Solver time in seconds.
  double time_io_ = 0.0;       ///< I/O time in seconds.

  // Initialization steps.
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

  // Conversion of option strings to enums.
  int GetPhysic(std::string physic_arg);
  utils::enums::implemType GetImplem(std::string implem_arg);
  utils::enums::methodType GetMethod(std::string method_arg);
  utils::enums::meshType GetMesh(std::string mesh_arg);
  model::AnisotropyType GetAnisotropy(std::string anisotropy_arg);
};

#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
