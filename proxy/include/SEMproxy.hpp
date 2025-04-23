//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include <SEMsolver.hpp>
#include <utils/argsparse.hpp>
#include <utils/utils.hpp>
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
  SEMproxy(int ex, int ey, int ez, float lx);

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

  // get information from mesh
  void getMeshInfo();

  SEMmesh myMesh;

  // Getter and setter for myRHSTerm
  arrayReal getMyRHSTerm() const;
  void setMyRHSTerm(const arrayReal &value);

  // Getter and setter for pnGlobal
  arrayReal getPnGlobal() const; 
  void setPnGlobal(const arrayReal &value); 

  // Getter and setter for rhsElement
  vectorInt getRhsElement() const;
  void setRhsElement(const vectorInt &value);

private:
  int i1 = 0;
  int i2 = 1;

  SEMinfo myInfo;

  SEMsolver mySolver;
  SolverUtils myUtils;

  // arrays
  arrayReal myRHSTerm;
  arrayReal pnGlobal;
  vectorInt rhsElement;

  // initialize source and RHS
  void init_source();

  // allocate arrays and vectors
  void init_arrays();
};

#endif /* SEMPROXY_HPP_ */
