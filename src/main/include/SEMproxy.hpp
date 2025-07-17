//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include <SEMsolver.hpp>
#include <argsparse.hpp>
#include <utils.hpp>
#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif

/**
 * @class SEMproxy
 */

class SEMproxy {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  SEMproxy(int argc, char *argv[]);

  SEMproxy(int ex, int ey, int ez, float lx)
      : myMesh(ex, ey, ez, lx, lx, lx, 2), mySolver(myMesh) {}

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

  void saveCtrlSlice(int iteration, int i);

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

private:
  int i1 = 0;
  int i2 = 1;

  const int myNumberOfRHS = 1;
  const float myTimeStep = 0.001;
  const float f0 = 10.;
  const float myTimeMax = 1.5;
  const int sourceOrder = 1;
  const static int order = 2;
  int myNumSamples = myTimeMax / myTimeStep;
  int myElementSource = 0;

  CartesianSEMmesh<float, int, int, order> myMesh;

  SEMsolver mySolver;
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
