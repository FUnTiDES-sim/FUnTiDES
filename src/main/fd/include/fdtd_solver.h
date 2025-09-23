//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef FDTD_SOLVER_HPP_
#define FDTD_SOLVER_HPP_

#include <memory>
#include <args_parse.h>
#include <utils.h>
#include "fdtd_grids.h"
#include "fdtd_stencils.h"
#include "fdtd_kernels.h"

/**
 * @class fdtd_solver
 */

class fdtd_solver {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  fdtd_solver() = default;
    
  fdtd_solver(fdtd_grids& grids, 
              fdtd_kernels& kernels, 
              fdtd_stencils& stencils): m_grids(grids), 
                                        m_kernels(kernels), 
                                        m_stencils(stencils){}

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~fdtd_solver(){};

  // compute one step
  void compute_one_step(int itime);

  
private:
  fdtd_options m_opt;

  int m_xsrc=100;
  int m_ysrc=100;
  int m_zsrc=100;

  fdtd_grids    & m_grids;
  fdtd_stencils & m_stencils;
  fdtd_kernels  & m_kernels;


};

#endif /* SEMPROXY_HPP_ */
