//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef FDPROXY_HPP_
#define FDPROXY_HPP_

#include <memory>
#include <args_parse.h>
#include <utils.h>
#include "FDproxyOptions.hpp"
#include "FDTDGrids.hpp"
#include "FDTDStencils.hpp"
#include "FDTDkernels.hpp"
#include "FDTDio.hpp"
#include "read_sepfile.h"

/**
 * @class FDproxy
 */

class FDproxy {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  FDproxy(const FDProxyOptions& opt);

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~FDproxy(){};

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void initFDElem();
  
  // initialize source and RHS
  void initSource();

  // compute one step
  void computeOneStep(int itime);

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

private:
  FDProxyOptions m_opt;
  
  int i1 = 0;
  int i2 = 1;

  int ncoefsX;
  int ncoefsY;
  int ncoefsZ;

  int   nSamples;
  float timeStep;
  float timeMax;
  
  int   sourceOrder;
  float f0;
  float vmin;
  float vmax;
  float lambdamax;

  int xsrc=-1;
  int ysrc=-1;
  int zsrc=-1;

  FDTDGRIDS    myGrids;
  FDTDStencils myStencils;
  FDTDKernels  myKernels;
  FDTDio       myIO;
  SolverUtils myUtils;

};

#endif /* SEMPROXY_HPP_ */
