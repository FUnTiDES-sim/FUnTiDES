//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include "SEMproxyOptions.hpp"
#include <solver_factory.h>
#include <argsparse.hpp>
#include <utils.hpp>
#include <memory>
#include <string>
#include <variant>
#include <model_struct.h>
#include <model_unstruct.h>

/**
 * @class SEMproxy
 */

class SEMproxy {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  SEMproxy(const SemProxyOptions& cfg);

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~SEMproxy(){};

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void initFiniteElem() {
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
  * Format: space-separated matrix with blank lines between rows for 3D plotting
  */
  void saveSlice(const VECTOR_REAL_VIEW& host_slice,
                int sizex, int sizey, const std::string& filepath) const;

  void saveSnapshot(int timesample) const;

  /**
  * @brief Computes optimal time step using CFL stability condition for seismic wave propagation
  *
  * Calculates the maximum stable time step for explicit finite difference schemes
  * using the CFL condition: dt ≤ CFL_factor × min_spacing / (√D × v_max)
  * where D is the number of dimension.
  *
  * @param cfl_factor Stability factor (0.5-0.7 for 2nd-order, 0.3-0.4 for higher-order schemes)
  * @return dt the max timestep.
  *
  * @note Typical values: 0.5-0.7 for 2nd-order, 0.3-0.4 for higher-order schemes
  * @warning Must be called before time-stepping loop to ensure numerical stability
  *
  */
  float find_cfl_dt(float cfl_factor);

private:
  int i1 = 0;
  int i2 = 1;

  // proper to cartesian mesh
  // or any structured mesh
  int nb_elements_[3] = {0};
  int nb_nodes_[3] = {0};

  // snapshots
  bool is_snapshots_;
  int snap_time_interval_;
  std::string snap_folder_;

  // time parameters
  float dt_;
  float timemax_;
  int num_sample_;
  // source parameters
  const int myNumberOfRHS = 1;
  const float f0 = 10.;
  const int sourceOrder = 2;
  int myElementSource = 0;

  std::variant <
    model::ModelStruct<float, int, 1>,
    model::ModelStruct<float, int, 2>,
    model::ModelStruct<float, int, 3>,
    model::ModelStruct<float, int, 4>,
    model::ModelUnstruct<float, int>
    > m_mesh_storage;
  model::ModelApi<float, int>* m_mesh;
  std::unique_ptr<SolverBase> m_solver;
  SolverUtils myUtils;

  // arrays
  arrayReal myRHSTerm;
  arrayReal pnGlobal;
  vectorInt rhsElement;
  arrayReal rhsWeights;

  // initialize source and RHS
  void init_source();

  // allocate arrays and vectors
  void init_arrays();


  // private methods to pars argv options
  int getPhysic ( string physicArg );
  SolverFactory::implemType getImplem ( string implemArg );
  SolverFactory::methodType getMethod ( string methodArg );
  SolverFactory::meshType getMesh ( string meshArg );
};

#endif /* SEMPROXY_HPP_ */
