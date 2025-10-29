//************************************************************************
// Finite Difference Time Domain (FDTD) Acoustic Simulation
// Version 0.0.1
//
// src/main/fd/fdtd_proxy.cpp
//
// Implementation of the FDTD proxy application main interface
//
// This file implements the FdtdProxy class, which orchestrates the entire
// FDTD simulation workflow including initialization, time-stepping, and
// output generation for acoustic wave propagation modeling.
//
// Copyright (c) 2025
// License: [Specify license here]
//************************************************************************

#include "fd_proxy.h"

#include <chrono>
#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

#include "data_type.h"

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::nanoseconds;
using std::chrono::system_clock;
using std::chrono::time_point;

namespace
{

/**
 * @brief Prints a section separator line.
 */
void PrintSeparator()
{
  std::cout << "--------------------------------------" << std::endl;
}

/**
 * @brief Converts microseconds to seconds.
 * @param time_us Time in microseconds
 * @return Time in seconds
 */
inline double MicrosecondsToSeconds(double time_us) { return time_us / 1E6; }

}  // namespace

namespace fdtd
{

FdtdProxy::FdtdProxy(const fdtd::options::FdtdOptions& opt)
    : opt_(opt),
      grids_(),
      stencils_(),
      kernels_(),
      io_(),
      utils_(),
      abckernels_(),
      solver_(grids_, kernels_, abckernels_, stencils_, source_receivers_)
{
}

void FdtdProxy::InitFdtd()
{
  std::cout << "+======================================" << std::endl;
  std::cout << "saveSnapshots=" << opt_.output.save_snapshots
            << " snapShotInterval=" << opt_.output.snapshot_interval
            << std::endl;
  PrintSeparator();
  std::cout << std::endl;

  InitializeGrid();
  InitializeStencils();
  InitializeVelocityModel();
  InitializeModelArrays();
  InitializeWavefieldArrays();
  InitializeSource();
  InitializeBoundaries();

  std::cout << "solver initialization done" << std::endl;
  PrintSeparator();
}

void FdtdProxy::InitializeGrid()
{
  std::cout << "geometry init" << std::endl;
  grids_.InitGrid(opt_);
  PrintSeparator();
  std::cout << "dx=" << grids_.dx() << " dy=" << grids_.dy()
            << " dz=" << grids_.dz() << std::endl;
  std::cout << "nx=" << grids_.nx() << " ny=" << grids_.ny()
            << " nz=" << grids_.nz() << std::endl;
}

void FdtdProxy::InitializeStencils()
{
  std::cout << "stencil init" << std::endl;
  PrintSeparator();
  stencils_.initStencilsCoefficients(opt_, grids_.dx(), grids_.dy(),
                                     grids_.dz());
  std::cout << "stencil coefficients" << std::endl;
  std::cout << "lx=" << stencils_.lx << " ly=" << stencils_.ly
            << " lz=" << stencils_.lz << std::endl;
  std::cout << "coef0=" << stencils_.coef0 << std::endl;

  for (int i = 0; i < stencils_.ncoefsX; i++)
  {
    std::cout << "coefx[" << i << "]=" << stencils_.coefx[i] << " ";
  }
  std::cout << std::endl;

  for (int i = 0; i < stencils_.ncoefsY; i++)
  {
    std::cout << "coefy[" << i << "]=" << stencils_.coefy[i] << " ";
  }
  std::cout << std::endl;

  for (int i = 0; i < stencils_.ncoefsZ; i++)
  {
    std::cout << "coefz[" << i << "]=" << stencils_.coefz[i] << " ";
  }
  std::cout << std::endl;
}

void FdtdProxy::InitializeVelocityModel()
{
  std::cout << std::endl;
  std::cout << "velocity model init" << std::endl;
  std::cout << "vmin=" << opt_.velocity.vmin << " vmax=" << opt_.velocity.vmax
            << std::endl;
  PrintSeparator();

  velocity_min_ = opt_.velocity.vmin;
  velocity_max_ = opt_.velocity.vmax;
  wavelength_max_ = opt_.velocity.vmax / (2.5 * opt_.source.f0);
  time_step_ = opt_.time.time_step;
  time_max_ = opt_.time.time_max;

  std::cout << "user defined time step=" << std::scientific << time_step_
            << std::endl;
  std::cout << "user defined max time=" << std::fixed << opt_.time.time_max
            << std::endl;
  PrintSeparator();

  // Compute time step from CFL condition if not user-defined
  if (time_step_ == 0.0f)
  {
    time_step_ = stencils_.compute_dt_sch(velocity_max_);
    std::cout << "compute time step from CFL condition" << std::endl;
  }
  else
  {
    std::cout << "user defined time step" << std::endl;
  }

  num_time_samples_ = static_cast<int>(time_max_ / time_step_);
  std::cout << "timeStep=" << std::scientific << time_step_ << std::endl;
  std::cout << "nSamples=" << num_time_samples_ << std::endl;
  PrintSeparator();
}

void FdtdProxy::InitializeModelArrays()
{
  std::cout << "model init" << std::endl;
  grids_.InitModelArrays(opt_);
  std::cout << "model init done" << std::endl;
  PrintSeparator();
}

void FdtdProxy::InitializeWavefieldArrays()
{
  kernels_.initFieldsArrays(grids_.nx(), grids_.ny(), grids_.nz(), stencils_.lx,
                            stencils_.ly, stencils_.lz);
  std::cout << "arrays init done" << std::endl;
  PrintSeparator();
}

void FdtdProxy::InitializeSource()
{
  source_frequency_ = opt_.source.f0;
  source_order_ = opt_.source.source_order;
  std::cout << "central freq and source order" << std::endl;
  std::cout << "f0=" << source_frequency_ << std::endl;
  std::cout << "sourceOrder=" << source_order_ << std::endl;
  PrintSeparator();

  // Set source position (use grid center if not specified)
  source_receivers_.xsrc =
      (opt_.source.xs < 0) ? grids_.nx() / 2 : opt_.source.xs;
  source_receivers_.ysrc =
      (opt_.source.ys < 0) ? grids_.ny() / 2 : opt_.source.ys;
  source_receivers_.zsrc =
      (opt_.source.zs < 0) ? grids_.nz() / 2 : opt_.source.zs;

  std::cout << "source position" << std::endl;
  std::cout << "xsrc=" << source_receivers_.xsrc
            << " ysrc=" << source_receivers_.ysrc
            << " zsrc=" << source_receivers_.zsrc << std::endl;
  PrintSeparator();

  InitSource();
  std::cout << "source init done" << std::endl;
  PrintSeparator();
}

void FdtdProxy::InitializeBoundaries()
{
  std::cout << "boundary init" << std::endl;
  // abckernels_.Initialize(opt_);
  if (opt_.boundary.use_sponge)
  {
    int nx = grids_.nx();
    int ny = grids_.ny();
    int nz = grids_.nz();
    abckernels_.spongeArray =
        allocateVector<vectorReal>(nx * ny * nz, "spongeArray");
    abckernels_.defineSpongeBoundary(nx, ny, nz);
    std::cout << "sponge boundary init done" << std::endl;
  }
  if (opt_.boundary.use_pml)
  {
    int nx = grids_.nx();
    int ny = grids_.ny();
    int nz = grids_.nz();
    int ndampx = grids_.ndampx();
    int ndampy = grids_.ndampy();
    int ndampz = grids_.ndampz();
    int x1 = grids_.x1();
    int x2 = grids_.x2();
    int x3 = grids_.x3();
    int x4 = grids_.x4();
    int x5 = grids_.x5();
    int x6 = grids_.x6();
    int y1 = grids_.y1();
    int y2 = grids_.y2();
    int y3 = grids_.y3();
    int y4 = grids_.y4();
    int y5 = grids_.y5();
    int y6 = grids_.y6();
    int z1 = grids_.z1();
    int z2 = grids_.z2();
    int z3 = grids_.z3();
    int z4 = grids_.z4();
    int z5 = grids_.z5();
    int z6 = grids_.z6();
    float dx = grids_.dx();
    float dy = grids_.dy();
    float dz = grids_.dz();
    float dt_sch = 0.001f;
    float vmax = opt_.velocity.vmax;

    printf("PML params: nx=%d ny=%d nz=%d ndampx=%d ndampy=%d ndampz=%d\n", nx,
           ny, nz, ndampx, ndampy, ndampz);
    // allocate eta array
    abckernels_.eta =
        allocateVector<vectorReal>((nx + 2) * (ny + 2) * (nz + 2), "eta");
    vectorReal& eta = abckernels_.eta;
    abckernels_.init_eta(nx, ny, nz, ndampx, ndampy, ndampz, x1, x2, x3, x4, x5,
                         x6, y1, y2, y3, y4, y5, y6, z1, z2, z3, z4, z5, z6, dx,
                         dy, dz, dt_sch, vmax, eta);
    std::cout << "PML boundary init done" << std::endl;
  }
  PrintSeparator();
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
  }
}

void FdtdProxy::Run()
{
  nanoseconds total_compute_time{0};
  nanoseconds total_output_time{0};

  for (int index_time_sample = 0; index_time_sample < num_time_samples_;
       index_time_sample++)
  {
    // Compute one time step
    auto start_compute_time = system_clock::now();
    if (opt_.boundary.use_sponge)
    {
      solver_.compute_one_stepSB(index_time_sample, time_index_current_,
                                 time_index_next_);
    }
    if (opt_.boundary.use_pml)
    {
      solver_.compute_one_stepPML(index_time_sample, time_index_current_,
                                  time_index_next_);
    }
    total_compute_time +=
        duration_cast<nanoseconds>(system_clock::now() - start_compute_time);

    // Output snapshots at specified intervals
    auto start_output_time = system_clock::now();
    if (index_time_sample % opt_.output.snapshot_interval == 0)
    {
      io_.outputPnValues(index_time_sample, time_index_current_, grids_,
                         kernels_, stencils_, opt_, source_receivers_);
    }

    // Swap time indices for next iteration
    std::swap(time_index_current_, time_index_next_);

    total_output_time +=
        duration_cast<nanoseconds>(system_clock::now() - start_output_time);
    std::cout.flush();
  }

  PrintPerformanceMetrics(total_compute_time, total_output_time);
}

void FdtdProxy::PrintPerformanceMetrics(
    const nanoseconds& total_compute_time,
    const nanoseconds& total_output_time) const
{
  const double kernel_time_us = static_cast<double>(
      duration_cast<microseconds>(total_compute_time).count());
  const double output_time_us = static_cast<double>(
      duration_cast<microseconds>(total_output_time).count());

  std::cout << "------------------------------------------------ " << std::endl;
  std::cout << "\n---- Elapsed Kernel Time : "
            << MicrosecondsToSeconds(kernel_time_us) << " seconds."
            << std::endl;
  std::cout << "---- Elapsed Output Time : "
            << MicrosecondsToSeconds(output_time_us) << " seconds."
            << std::endl;
  std::cout << "------------------------------------------------ " << std::endl;
}

}  // namespace fdtd
