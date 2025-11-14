//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef FDTD_SOLVER_HPP_
#define FDTD_SOLVER_HPP_

#include <fd_grids.h>
#include <utils.h>

#include <memory>

#include "fd_abckernels.h"
#include "fd_kernels.h"
#include "fd_source_receivers.h"
#include "fd_stencils.h"

/**
 * @class fdtd_solver
 */

class FdtdSolver
{
 public:
  /**
   * @brief Constructor of the SEMproxy class
   */

  FdtdSolver(model::fdgrid::FdtdGrids& grids,
             fdtd::kernel::FdtdKernels& kernels,
             fdtd::abckernel::FdtdAbcKernels& abckernels,
             fdtd::stencils::FdtdStencils& stencils,
             FdtdSourceReceivers& source_receivers)
      : m_grids(grids),
        m_kernels(kernels),
        m_abckernels(abckernels),
        m_stencils(stencils),
        m_source_receivers(source_receivers){};

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~FdtdSolver(){};

  // compute one step with sponge boundary
  void compute_one_stepSB(int itime, int i1, int i2);

  // compute one step with PML
  void compute_one_stepPML(int itime, int i1, int i2);

 private:
  model::fdgrid::FdtdGrids& m_grids;
  fdtd::kernel::FdtdKernels& m_kernels;
  fdtd::abckernel::FdtdAbcKernels& m_abckernels;
  fdtd::stencils::FdtdStencils& m_stencils;
  FdtdSourceReceivers& m_source_receivers;
};

#endif /* FDTD_SOLVER_HPP_ */
