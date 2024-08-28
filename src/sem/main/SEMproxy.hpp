//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include "SEMsolver.hpp"
#include "utils.hpp"

/**
 * @class SEMproxy
 */

class SEMproxy
{
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
  void initFiniteElem();

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

  // initialize source and RHS
  void init_source();

  // allocate arrays and vectors
  void init_arrays();

  // get information from mesh
  void getMeshInfo();


private:

  int i1=0;
  int i2=1;

  SEMinfo myInfo;
  SEMmesh myMesh;

  SEMsolver mySolver;
  SolverUtils myUtils;

  // arrays
  arrayReal myRHSTerm;
  arrayReal pnGlobal;
  vectorInt rhsElement;

};

#endif /* SEMPROXY_HPP_ */
