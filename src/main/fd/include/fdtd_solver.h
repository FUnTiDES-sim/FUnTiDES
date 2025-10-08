//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef FDTD_SOLVER_HPP_
#define FDTD_SOLVER_HPP_

#include <args_parse.h>
#include <fdtd_grids.h>
#include <utils.h>

#include <memory>

#include "fdtd_kernels.h"
#include "fdtd_source_receivers.h"
#include "fdtd_stencils.h"

/**
 * @class fdtd_solver
 */

class FdtdSolver
{
 public:
  /**
   * @brief Constructor of the SEMproxy class
   */

  FdtdSolver(model::fdgrid::FdtdGrids& grids, FdtdKernels& kernels,
             FdtdStencils& stencils, FdtdSourceReceivers& source_receivers)
      : m_grids(grids),
        m_kernels(kernels),
        m_stencils(stencils),
        m_source_receivers(source_receivers)
  {
  }

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~FdtdSolver(){};

  // compute one step
  void compute_one_step(int itime, int i1, int i2);

 private:
  model::fdgrid::FdtdGrids& m_grids;
  FdtdStencils& m_stencils;
  FdtdKernels& m_kernels;
  FdtdSourceReceivers& m_source_receivers;
};

#endif /* FDTD_SOLVER_HPP_ */
