
//************************************************************************
// Finite Difference Time Domain (FDTD) Acoustic Simulation
//
// src/main/fd/fdtd_proxy.cpp
//
// Implementation of the FDTD proxy application main interface
//
// This file implements the FdtdProxy class, which orchestrates the entire
// FDTD simulation workflow including initialization, time-stepping, and
// output generation for acoustic wave propagation modeling.
//************************************************************************

#include "fdtd_proxy.h"

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

#include "data_type.h"

FdtdProxy::FdtdProxy(const FdtdOptions& opt)
    : opt_(opt),
      grids_(),
      stencils_(),
      kernels_(),
      io_(),
      utils_(),
      solver_(grids_, kernels_, stencils_, source_receivers_)
{
}

void FdtdProxy::InitFdtd()
{
  printf("+======================================\n");
  printf("saveSnapshots=%d snapShotInterval=%d\n", opt_.output.save_snapshots,
         opt_.output.snapshot_interval);
  printf("--------------------------------------\n");
  printf("\n");

  // Initialize computational grid geometry
  printf("geometry init\n");
  grids_.InitGrid(opt_);
  printf("--------------------------------------\n");
  printf("dx=%f dy=%f dz=%f\n", grids_.dx(), grids_.dy(), grids_.dz());
  printf("nx=%d ny=%d nz=%d\n", grids_.nx(), grids_.ny(), grids_.nz());

  // Initialize finite difference stencil coefficients
  printf("stencil init\n");
  printf("--------------------------------------\n");
  stencils_.initStencilsCoefficients(opt_, grids_.dx(), grids_.dy(),
                                     grids_.dz());
  printf("stencil coefficients\n");
  printf("lx=%d ly=%d lz=%d\n", stencils_.lx, stencils_.ly, stencils_.lz);
  printf("coef0=%f\n", stencils_.coef0);
  for (int i = 0; i < stencils_.ncoefsX; i++)
  {
    printf("coefx[%d]=%f ", i, stencils_.coefx[i]);
  }
  printf("\n");
  for (int i = 0; i < stencils_.ncoefsY; i++)
  {
    printf("coefy[%d]=%f ", i, stencils_.coefy[i]);
  }
  printf("\n");
  for (int i = 0; i < stencils_.ncoefsZ; i++)
  {
    printf("coefz[%d]=%f ", i, stencils_.coefz[i]);
  }
  printf("\n");

  // Initialize velocity model parameters
  printf("\n");
  printf("velocity model init\n");
  printf("vmin=%f vmax=%f\n", opt_.velocity.vmin, opt_.velocity.vmax);
  printf("--------------------------------------\n");
  velocity_min_ = opt_.velocity.vmin;
  velocity_max_ = opt_.velocity.vmax;
  wavelength_max_ = opt_.velocity.vmax / (2.5 * opt_.source.f0);
  time_step_ = opt_.time.time_step;
  time_max_ = opt_.time.time_max;
  printf("user defined time step=%e\n", time_step_);
  printf("user defined max time=%f\n", opt_.time.time_max);
  printf("--------------------------------------\n");

  // Compute time step from CFL condition if not user-defined
  if (time_step_ == 0)
  {
    time_step_ = stencils_.compute_dt_sch(velocity_max_);
    printf("compute time step from CFL condition\n");
  }
  else
  {
    printf("user defined time step\n");
  }
  num_time_samples_ = time_max_ / time_step_;
  printf("timeStep=%e\n", time_step_);
  printf("nSamples=%d\n", num_time_samples_);
  printf("--------------------------------------\n");

  // Initialize model arrays (velocity, density, etc.)
  printf("model init\n");
  grids_.InitModelArrays(opt_);
  printf("model init done\n");
  printf("--------------------------------------\n");

  // Allocate and initialize wavefield arrays
  kernels_.initFieldsArrays(grids_.nx(), grids_.ny(), grids_.nz(), stencils_.lx,
                            stencils_.ly, stencils_.lz);
  printf("arrays init done\n");
  printf("--------------------------------------\n");

  // Configure seismic source parameters
  source_frequency_ = opt_.source.f0;
  source_order_ = opt_.source.source_order;
  printf("central freq and source order\n");
  printf("f0=%f\n", source_frequency_);
  printf("sourceOrder=%d\n", source_order_);
  printf("--------------------------------------\n");

  // Set source position (use grid center if not specified)
  source_receivers_.xsrc = opt_.source.xs;
  source_receivers_.ysrc = opt_.source.ys;
  source_receivers_.zsrc = opt_.source.zs;
  if (source_receivers_.xsrc < 0)
  {
    source_receivers_.xsrc = grids_.nx() / 2;
  }
  if (source_receivers_.ysrc < 0)
  {
    source_receivers_.ysrc = grids_.ny() / 2;
  }
  if (source_receivers_.zsrc < 0)
  {
    source_receivers_.zsrc = grids_.nz() / 2;
  }
  printf("source position\n");
  printf("xsrc=%d ysrc=%d zsrc=%d\n", source_receivers_.xsrc,
         source_receivers_.ysrc, source_receivers_.zsrc);
  printf("--------------------------------------\n");

  InitSource();
  printf("source init done\n");
  printf("--------------------------------------\n");

  // Define absorbing boundary conditions (sponge layers)
  kernels_.defineSpongeBoundary(grids_.nx(), grids_.ny(), grids_.nz());
  printf("sponge boundary init done\n");
  printf("--------------------------------------\n");

  printf("solver initialization done\n");
  printf("--------------------------------------\n");
}

void FdtdProxy::InitSource()
{
  // Compute source term (e.g., Ricker wavelet)
  kernels_.RHSTerm = allocateVector<vectorReal>(num_time_samples_, "RHSTerm");

  std::vector<float> source_term = utils_.computeSourceTerm(
      num_time_samples_, time_step_, source_frequency_, source_order_);
  for (int i = 0; i < num_time_samples_; i++)
  {
    kernels_.RHSTerm[i] = source_term[i];
    // std::cout << "sample " << i << "\t: sourceTerm = " << source_term[i]
    //           << std::endl;
  }
}

void FdtdProxy::Run()
{
  time_point<system_clock> start_compute_time, start_output_time,
      total_compute_time, total_output_time;

  for (int index_time_sample = 0; index_time_sample < num_time_samples_;
       index_time_sample++)
  {
    // Compute one time step
    start_compute_time = system_clock::now();
    solver_.compute_one_step(index_time_sample, time_index_current_,
                             time_index_next_);
    total_compute_time += system_clock::now() - start_compute_time;

    // Output snapshots at specified intervals
    start_output_time = system_clock::now();
    if (index_time_sample % opt_.output.snapshot_interval == 0)
    {
      io_.outputPnValues(index_time_sample, time_index_current_, grids_,
                         kernels_, stencils_, opt_, source_receivers_);
    }

    // Swap time indices for next iteration
    std::swap(time_index_current_, time_index_next_);

    total_output_time += system_clock::now() - start_output_time;
    fflush(stdout);
  }

  // Report performance metrics
  float kernel_time_ms = time_point_cast<microseconds>(total_compute_time)
                             .time_since_epoch()
                             .count();
  float output_time_ms = time_point_cast<microseconds>(total_output_time)
                             .time_since_epoch()
                             .count();

  std::cout << "------------------------------------------------ " << std::endl;
  std::cout << "\n---- Elapsed Kernel Time : " << kernel_time_ms / 1E6
            << " seconds." << std::endl;
  std::cout << "---- Elapsed Output Time : " << output_time_ms / 1E6
            << " seconds." << std::endl;
  std::cout << "------------------------------------------------ " << std::endl;
}
