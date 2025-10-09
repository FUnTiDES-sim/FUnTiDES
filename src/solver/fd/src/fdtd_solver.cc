//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "fdtd_solver.h"

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

#include "data_type.h"

// compute one time step
void FdtdSolver::compute_one_step(int itime, int i1, int i2)
{
  int x3 = m_stencils.lx;
  int x4 = m_grids.nx() - m_stencils.lx;
  int y3 = m_stencils.ly;
  int y4 = m_grids.ny() - m_stencils.ly;
  int z3 = m_stencils.lz;
  int z4 = m_grids.nz() - m_stencils.lz;
  m_kernels.addRHS(itime, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                   m_stencils.lx, m_stencils.ly, m_stencils.lz,
                   m_source_receivers.xsrc, m_source_receivers.ysrc,
                   m_source_receivers.zsrc, m_grids.vp(), m_kernels.RHSTerm,
                   m_kernels.pnGlobal);

  // printf("addRHS done\n");
  FDFENCE
  // inner points

  m_kernels.inner3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                    m_stencils.lx, m_stencils.ly, m_stencils.lz, x3, x4, y3, y4,
                    z3, z4, m_stencils.coef0, m_stencils.coefx,
                    m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                    m_kernels.pnGlobal);
  // printf("inner3D done\n");
  FDFENCE
  // apply sponge boundary to wavefield
  m_kernels.applySponge(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                        m_stencils.lx, m_stencils.ly, m_stencils.lz, x3, x4, y3,
                        y4, z3, z4, m_kernels.spongeArray, m_kernels.pnGlobal);
  // printf("applySponge done\n");
  FDFENCE
}

//************************************************************************
// End of file
//************************************************************************
