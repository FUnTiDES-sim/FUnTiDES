//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef FDTD_PROXY_HPP_
#define FDTD_PROXY_HPP_

#include <memory>
#include <args_parse.h>
#include <utils.h>
#include "fdtd_solver.h"
#include "fdtd_options.h"
#include "fdtd_grids.h"
#include "fdtd_stencils.h"
#include "fdtd_kernels.h"
#include "fdtd_io.h"
#include "read_sepfile.h"

/**
 * @class FDproxy
 */

class fdtd_proxy {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  fdtd_proxy(const fdtd_options & opt);

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~fdtd_proxy(){};

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void init_fdtd(); 

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

private:
  fdtd_options m_opt;
  
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

  // initialize source and RHS
  void init_source();

  fdtd_grids    m_grids;
  fdtd_stencils m_stencils;
  fdtd_kernels  m_kernels;
  fdtd_io       m_io;
  SolverUtils   m_utils;
  fdtd_solver   m_solver;

};

#endif /* SEMPROXY_HPP_ */
