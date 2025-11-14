#include "fd_solver.h"

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

#include "data_type.h"

/**
 * @brief Computes a single time step in the FDTD acoustic wave simulation.
 *
 * This function performs one complete iteration of the Finite-Difference
 * Time-Domain (FDTD) algorithm for acoustic wave propagation. It executes
 * three main operations in sequence: adding source terms, computing inner
 * domain updates, and applying absorbing boundary conditions.
 *
 * @param itime Current time step index in the simulation
 * @param i1 First buffer index for the two-buffer time stepping scheme
 * @param i2 Second buffer index for the two-buffer time stepping scheme
 *
 * @details The computation proceeds in three phases:
 *
 * 1. Source Term Addition: Injects acoustic energy from sources into
 *    the pressure field using the RHS (right-hand side) term at specified
 *    source locations.
 *
 * 2. Inner Domain Update: Applies the finite-difference stencil to
 *    compute spatial derivatives and update the pressure field in the
 *    interior domain, excluding the boundary layers defined by the stencil
 *    half-lengths (lx, ly, lz).
 *
 * 3. Sponge Boundary Application: Attenuates outgoing waves at the
 *    domain boundaries to prevent non-physical reflections using absorbing
 *    boundary conditions.
 *
 * @note The function uses FDFENCE macros between operations to ensure proper
 *       memory synchronization in parallel execution environments.
 *
 * @warning This function assumes all grid arrays and stencil coefficients
 *          have been properly initialized before invocation.
 *
 * @see FdtdSolver::addRHS()
 * @see FdtdSolver::inner3D()
 * @see FdtdSolver::applySponge()
 */
void FdtdSolver::compute_one_stepSB(int itime, int i1, int i2)
{
  m_kernels.addRHS(itime, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                   m_stencils.lx, m_stencils.ly, m_stencils.lz,
                   m_source_receivers.xsrc, m_source_receivers.ysrc,
                   m_source_receivers.zsrc, m_grids.vp(), m_kernels.RHSTerm,
                   m_kernels.pnGlobal);

  // printf("addRHS done\n");
  FDFENCE
  // inner points

  m_kernels.inner3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                    m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x3(),
                    m_grids.x4(), m_grids.y3(), m_grids.y4(), m_grids.z3(),
                    m_grids.z4(), m_stencils.coef0, m_stencils.coefx,
                    m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                    m_kernels.pnGlobal);
  // printf("inner3D done\n");
  FDFENCE
  // apply sponge boundary to wavefield
  m_kernels.applySponge(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                        m_stencils.lx, m_stencils.ly, m_stencils.lz,
                        m_grids.x3(), m_grids.x4(), m_grids.y3(), m_grids.y4(),
                        m_grids.z3(), m_grids.z4(), m_abckernels.spongeArray,
                        m_kernels.pnGlobal);
  // printf("applySponge done\n");
  FDFENCE
}

void FdtdSolver::compute_one_stepPML(int itime, int i1, int i2)
{
  m_kernels.addRHS(itime, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                   m_stencils.lx, m_stencils.ly, m_stencils.lz,
                   m_source_receivers.xsrc, m_source_receivers.ysrc,
                   m_source_receivers.zsrc, m_grids.vp(), m_kernels.RHSTerm,
                   m_kernels.pnGlobal);

  // printf("addRHS done\n");
  FDFENCE

  // update PML
  // up
  m_kernels.pml3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                  m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x1(),
                  m_grids.x6(), m_grids.y1(), m_grids.y6(), m_grids.z1(),
                  m_grids.z2(), m_stencils.coef0, m_grids.hdx_2(),
                  m_grids.hdy_2(), m_grids.hdz_2(), m_stencils.coefx,
                  m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                  m_abckernels.eta, m_kernels.phi, m_kernels.pnGlobal);
  // update front
  m_kernels.pml3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                  m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x1(),
                  m_grids.x6(), m_grids.y1(), m_grids.y2(), m_grids.z3(),
                  m_grids.z4(), m_stencils.coef0, m_grids.hdx_2(),
                  m_grids.hdy_2(), m_grids.hdz_2(), m_stencils.coefx,
                  m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                  m_abckernels.eta, m_kernels.phi, m_kernels.pnGlobal);
  // update left
  m_kernels.pml3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                  m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x1(),
                  m_grids.x2(), m_grids.y3(), m_grids.y4(), m_grids.z3(),
                  m_grids.z4(), m_stencils.coef0, m_grids.hdx_2(),
                  m_grids.hdy_2(), m_grids.hdz_2(), m_stencils.coefx,
                  m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                  m_abckernels.eta, m_kernels.phi, m_kernels.pnGlobal);
  // inner points
  m_kernels.inner3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                    m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x3(),
                    m_grids.x4(), m_grids.y3(), m_grids.y4(), m_grids.z3(),
                    m_grids.z4(), m_stencils.coef0, m_stencils.coefx,
                    m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                    m_kernels.pnGlobal);
  // update right
  m_kernels.pml3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                  m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x5(),
                  m_grids.x6(), m_grids.y3(), m_grids.y4(), m_grids.z3(),
                  m_grids.z4(), m_stencils.coef0, m_grids.hdx_2(),
                  m_grids.hdy_2(), m_grids.hdz_2(), m_stencils.coefx,
                  m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                  m_abckernels.eta, m_kernels.phi, m_kernels.pnGlobal);
  // update back
  m_kernels.pml3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                  m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x1(),
                  m_grids.x6(), m_grids.y5(), m_grids.y6(), m_grids.z3(),
                  m_grids.z4(), m_stencils.coef0, m_grids.hdx_2(),
                  m_grids.hdy_2(), m_grids.hdz_2(), m_stencils.coefx,
                  m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                  m_abckernels.eta, m_kernels.phi, m_kernels.pnGlobal);
  // update bottom
  m_kernels.pml3D(i1, i2, m_grids.nx(), m_grids.ny(), m_grids.nz(),
                  m_stencils.lx, m_stencils.ly, m_stencils.lz, m_grids.x1(),
                  m_grids.x6(), m_grids.y1(), m_grids.y6(), m_grids.z5(),
                  m_grids.z6(), m_stencils.coef0, m_grids.hdx_2(),
                  m_grids.hdy_2(), m_grids.hdz_2(), m_stencils.coefx,
                  m_stencils.coefy, m_stencils.coefz, m_grids.vp(),
                  m_abckernels.eta, m_kernels.phi, m_kernels.pnGlobal);
  FDFENCE
}
