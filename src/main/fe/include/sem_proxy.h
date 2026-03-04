#ifndef FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
#define FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
#include <data_type.h>
#include <utils.h>

#include <memory>
#include <string>
#include <variant>

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

using namespace solver::fe::enums;

/**
 * @class SEMproxy
 */
class SEMproxy
{
 public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  SEMproxy(const SemProxyOptions& cfg);

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~SEMproxy() {}

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void initFiniteElem()
  {
    init_arrays();
    init_source();
  };

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

  /**
   * Save slice in gnuplot matrix format (default - best for gnuplot)
   * Format: space-separated matrix with blank lines between rows for 3D
   * plotting
   */
  void saveSlice(const VECTOR_REAL_VIEW& host_slice, int sizex, int sizey,
                 const std::string& filepath) const;

  void saveSnapshot(int timesample, VECTOR_REAL_VIEW data) const;

  /**
   * @brief Computes optimal time step using CFL stability condition for seismic
   * wave propagation
   *
   * Calculates the maximum stable time step for explicit finite difference
   * schemes using the CFL condition: dt ≤ CFL_factor × min_spacing / (√D ×
   * v_max) where D is the number of dimension.
   *
   * @param cfl_factor Stability factor (0.5-0.7 for 2nd-order, 0.3-0.4 for
   * higher-order schemes)
   * @return dt the max timestep.
   *
   * @note Typical values: 0.5-0.7 for 2nd-order, 0.3-0.4 for higher-order
   * schemes
   * @warning Must be called before time-stepping loop to ensure numerical
   * stability
   *
   */
  float find_cfl_dt(float cfl_factor);

  void saveSnapshot(int timestep);

 private:
  // Domain Decomposition Parameter
  // Mocking MPI rank and size for now.
  // In a real MPI application, these would come from MPI_Comm_rank/size.
  model::CartesianParams<float, int> m_localParams;
  utils::DistributedContext dist_ctx_;
  utils::ParallelTopology par_topology_;

  // proper to cartesian mesh
  // or any structured mesh
  int nb_elements_[3] = {0};
  int nb_nodes_[3] = {0};
  std::array<float, 3> src_coord_ = {0};
  std::array<float, 3> rcv_coord_ = {0};
  float domain_size_[3] = {0};

  // snapshots
  bool is_snapshots_;
  int snap_time_interval_;
  std::string snap_folder_;

  // physics
  bool isElastic_;
  bool isAcoustoElastic_;
  bool freeSurface_;

  // time parameters
  float dt_;
  float timemax_;
  int num_sample_;
  // source parameters
  const int myNumberOfRHS = 1;
  const float f0 = 5.0f;
  const int sourceOrder = 2;
  int myElementSource = 0;

  std::shared_ptr<model::ModelApi<float, int>> m_mesh;
  std::unique_ptr<solver::fe::Solver> m_solver;

  // Boundary Synchronizer for DD
  std::unique_ptr<solver::fe::BoundarySynchronizer> m_syncer;

  SolverUtils myUtils;

  // arrays
  arrayReal myRHSTerm;
  vectorReal pnGlobalPrev;
  vectorReal pnGlobalCurr;
  vectorInt rhsElement;
  vectorInt rhsElementRcv;
  arrayReal rhsWeights;
  arrayReal rhsWeightsRcv;
  arrayReal pnAtReceiver;

  // elastic arrays
  arrayReal myRHSTermx;
  arrayReal myRHSTermy;
  arrayReal myRHSTermz;
  vectorReal uxnGlobalPrev;
  vectorReal uynGlobalPrev;
  vectorReal uznGlobalPrev;
  vectorReal uxnGlobalCurr;
  vectorReal uynGlobalCurr;
  vectorReal uznGlobalCurr;
  arrayReal uxnAtReceiver;
  arrayReal uynAtReceiver;
  arrayReal uznAtReceiver;

  // io controller
  std::shared_ptr<SemIOController> io_ctrl_;

  // initialize source and RHS
  void init_source();

  // allocate arrays and vectors
  void init_arrays();

  // private methods to pars argv options
  int getPhysic(string physicArg);
  implemType getImplem(string implemArg);
  methodType getMethod(string methodArg);
  meshType getMesh(string meshArg);
  model::AnisotropyType getAnisotropy(std::string anisotropyArg);
};
#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_H_
