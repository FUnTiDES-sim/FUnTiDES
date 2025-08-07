//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include "SolverBase.hpp"
#include <argsparse.hpp>
#include <utils.hpp>
#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif

#include <memory>

/**
 * @class SEMproxy
 */

class SEMproxy {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  SEMproxy(int argc, char *argv[]);

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~SEMproxy(){};

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void initFiniteElem() {
    // allocate arrays and vectors
    init_arrays();
    // initialize source and RHS
    init_source();
  };

  // void saveCtrlSlice(int iteration, int i);

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
                int size, const std::string& filepath);

private:
  int i1 = 0;
  int i2 = 1;

  // proper to cartesian mesh
  // or any structured mesh
  int nb_elements[3] = {0};
  int nb_nodes[3] = {0};

  Mesh myMesh;

  float m_dt;
  float m_maxTime;
  int m_numSamples;

  float m_f0;
  int m_sourceOrder;
  int m_elementSource;
  int m_numberOfRHS;

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
};

#endif /* SEMPROXY_HPP_ */
