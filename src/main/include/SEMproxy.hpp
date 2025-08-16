//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include "SEMproxyOptions.hpp"
#include "solverFactory.hpp"
#include <argsparse.hpp>
#include <utils.hpp>
#include <memory>

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

  const int myNumberOfRHS = 1;
  const float myTimeStep = 0.001;
  const float f0 = 10.;
  const float myTimeMax = 1.5;
  const int sourceOrder = 2;

  int myNumSamples = myTimeMax / myTimeStep;
  int myElementSource = 0;

  BaseMesh<float, int> const* myMesh = nullptr;

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
