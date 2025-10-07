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
#include <fdtd_grids.h>
#include "fdtd_stencils.h"
#include "fdtd_kernels.h"
#include "fdtd_source_receivers.h"

/**
 * @class fdtd_solver
 */

class fdtd_solver {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
    
  fdtd_solver(model::fdgrid::FdtdGrids& grids,
              fdtd_kernels& kernels, 
              fdtd_stencils& stencils,
              fdtd_source_receivers &source_receivers): m_grids(grids), 
                                        m_kernels(kernels), 
                                        m_stencils(stencils),
                                        m_source_receivers(source_receivers){}

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~fdtd_solver(){};

  // compute one step
  void compute_one_step(int itime,int i1,int i2);

  
private:
  model::fdgrid::FdtdGrids    & m_grids;
  fdtd_stencils & m_stencils;
  fdtd_kernels  & m_kernels;
  fdtd_source_receivers & m_source_receivers;


};

#endif /* FDTD_SOLVER_HPP_ */
