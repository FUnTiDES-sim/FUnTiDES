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

  float m_dt;
  float m_maxTime;
  int m_numSamples;

  float m_f0;
  int m_sourceOrder;
  int m_elementSource;
  int m_numberOfRHS;

  int m_numberOfNodes;
  int m_numberOfElements;
  int m_numberOfPointsPerElement;
  int m_numberOfInteriorNodes;

  std::unique_ptr<SolverBase> m_solver;
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
