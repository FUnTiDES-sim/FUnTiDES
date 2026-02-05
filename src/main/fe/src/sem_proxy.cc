//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <boundary_synchronizer.h>
#include <cartesian_partitioner.h>
#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <variant>

#ifdef USE_MPI
#include "mpi_backend.h"
#endif
#include "sem_solver.h"
#include "topology_factory.h"

using namespace SourceAndReceiverUtils;
using namespace solver::fe;
using namespace solver::fe::enums;
using namespace solver::fe::solver_factory;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  // Init MPI parameters
  int mpi_initialized = 0;
  init_mpi(&mpi_initialized);

  init_sim_params(opt);
  init_mesh_params(opt);
  init_topology();
  init_sync();
  init_time_params(opt);

  const methodType methodType = getMethod(opt.method);
  const implemType implemType = getImplem(opt.implem);
  const meshType meshType = getMesh(opt.mesh);
  const modelLocationType modelLocation = opt.isModelOnNodes
                                              ? modelLocationType::kOnNodes
                                              : modelLocationType::kOnElements;
  const physicType physicType =
      opt.isElastic ? physicType::kElastic : physicType::kAcoustic;

  m_solver = createSolver(methodType, implemType, meshType, modelLocation,
                          physicType, opt.order);

  // Setup Sponge Parameters
  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};

  // Note: m_solver->computeFEInit is now called in run() to pass partition
  // info. We manually call init arrays here if needed, but computeFEInit does
  // it too. For consistency with old code structure, we prep arrays now.
  initFiniteElem();

  // Prepare ADIOS2 dimensions
  // Global Dimensions
  // Note: For unstructured/general meshes, this is complex.
  // For Cartesian Struct Mesh:
  // Global Size = Global_Elements * Order + 1
  // Local Size = Local_Elements * Order + 1
  // Offset = Origin_Index (calculated from origin_x)
  size_t global_nx =
      m_localParams.global_lx / (m_localParams.global_lx / opt.ex) * opt.order +
      1;  // Approximation if uniform
  // More precise:
  size_t g_ex = static_cast<size_t>(opt.ex);
  size_t g_ey = static_cast<size_t>(opt.ey);
  size_t g_ez = static_cast<size_t>(opt.ez);

  // Just 1D decomposition on X supported by partitioner currently
  size_t g_nx = g_ex * opt.order + 1;
  size_t g_ny = g_ey * opt.order + 1;
  size_t g_nz = g_ez * opt.order + 1;
  std::vector<size_t> global_dims = {g_nx, g_ny, g_nz};

  // Calculate offset based on origin
  // origin_x is physical coordinate. We need node index.
  // dx = lx / ex.
  float dx = m_localParams.global_lx / opt.ex;
  size_t elem_offset_x =
      static_cast<size_t>(std::round(m_localParams.origin_x / dx));
  size_t node_offset_x = elem_offset_x * opt.order;

  std::vector<size_t> start_offsets = {node_offset_x, 0, 0};

  // Local Size
  // Note: ADIOS2 overlap handling.
  // Nodes are shared. If we write all local nodes, boundary nodes are written
  // twice. ADIOS2 allows this, but we must ensure we are consistent. Simple
  // approach: Write full local block.
  size_t l_nx = static_cast<size_t>(nb_nodes_[0]);
  size_t l_ny = static_cast<size_t>(nb_nodes_[1]);
  size_t l_nz = static_cast<size_t>(nb_nodes_[2]);
  std::vector<size_t> local_dims = {l_nx, l_ny, l_nz};

  // Refine Global Dims vector
  global_dims = {g_nx, g_ny, g_nz};

  io_ctrl_ = std::make_shared<SemIOController>(
      global_dims, start_offsets, local_dims, static_cast<size_t>(num_sample_),
      static_cast<size_t>(1));

  // snapshots settings
  is_snapshots_ = opt.snapshots;
  if (is_snapshots_)
  {
    snap_time_interval_ = opt.snap_time_interval;
  }
}

void SEMproxy::run()
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  bool isElastic = isElastic_;

  // Sponge params from options
  const float spongex =
      0;  // Configured earlier but variable scope issue in original
  const std::array<float, 3> sponge_size = {0, 0, 0};
  const bool surface_sponge = false;
  const float taper_delta = 0.015;

  // Initialize Solver with Partition Info & Compute Local Mass
  m_solver->computeFEInit(*m_mesh, sponge_size, surface_sponge, taper_delta);

  // Synchronize Mass Matrix (Critical for DD)
  if (par_topology_.isDistributed())
  {
    m_syncer->synchronize(m_solver->getMassMatrix(), par_topology_);
  }

  auto& M = m_solver->getMassMatrix();
  // Get the global node index of the first node of the source element
  int debugNodeIdx = m_mesh->globalNodeIndex(myElementSource, 0, 0, 0);

  if (!isElastic)
  {
    WavefieldAcoustic wavefield(pnGlobalPrev, pnGlobalCurr);
    RhsAcoustic rhs(myRHSTerm, rhsElement, rhsWeights);
    SEMsolverDataAcoustic solverData(wavefield, rhs);

    for (int indexTimeSample = 0; indexTimeSample < num_sample_;
         indexTimeSample++)
    {
      startComputeTime = system_clock::now();

      // Compute Local Forces
      m_solver->computeForces(dt_, indexTimeSample, solverData);

      // Synchronize Forces
      if (par_topology_.isDistributed())
      {
        for (int c = 0; c < m_solver->getNumComponents(); ++c)
        {
          m_syncer->synchronize(m_solver->getForceVector(c), par_topology_);
        }
      }

      m_solver->updateSolution(dt_, solverData);

      totalComputeTime += system_clock::now() - startComputeTime;
      startOutputTime = system_clock::now();

      if (indexTimeSample % 50 == 0)
      {
        m_solver->outputSolutionValues(indexTimeSample, rhsElement[0],
                                       pnGlobalPrev, "pnGlobal");
      }

      // Save slice in dat format
      if (is_snapshots_ && indexTimeSample % snap_time_interval_ == 0)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "Save snapshot." << std::endl;
        saveSnapshot(indexTimeSample, pnGlobalPrev);
      }

      // Save pressure at receiver
      const int order = m_mesh->getOrder();

      float varnp1 = 0.0;
      for (int i = 0; i < order + 1; i++)
      {
        for (int j = 0; j < order + 1; j++)
        {
          for (int k = 0; k < order + 1; k++)
          {
            int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
            int globalNodeOnElement =
                i + j * (order + 1) + k * (order + 1) * (order + 1);
#!/usr/bin/env python3
"""
Simple BP5 Visualization Script
Reads BP5 directory and plots all timesteps
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import argparse


def read_bp5_file(bp_dir, nx, ny, nz):
    """Read all timesteps from BP5 data.0 file"""

    data_file = os.path.join(bp_dir, "data.0")

    if not os.path.exists(data_file):
        print(f"Error: {data_file} not found")
        return None, None

    print(f"Reading {data_file}...")

    # Read raw binary file
    with open(data_file, 'rb') as f:
        raw_data = f.read()

    # Convert to float32 array
    num_floats = len(raw_data) // 4
    data_array = np.frombuffer(raw_data, dtype=np.float32, count=num_floats)

    print(f"Total floats read: {num_floats}")
    print(f"Grid size: {nx} x {ny} x {nz} = {nx*ny*nz} per timestep")

    # Calculate timesteps
    values_per_step = nx * ny * nz
    num_steps = num_floats // values_per_step
    remainder = num_floats % values_per_step

    print(f"Timesteps: {num_steps} complete + {remainder} remainder")

    # Extract all timesteps
    data_list = []
    for step in range(num_steps):
        start = step * values_per_step
        end = start + values_per_step
        step_data = data_array[start:end].reshape((nx, ny, nz))
        data_list.append(step_data)

    print(f"Extracted {len(data_list)} timesteps")
    return data_list, np.arange(len(data_list))


def plot_slices(data_3d, timestep, output_file, vmin=None, vmax=None):
    """Plot XY, XZ, YZ slices"""

    nx, ny, nz = data_3d.shape
    mid_x = nx // 2
    mid_y = ny // 2
    mid_z = nz // 2

    # Get slices
    slice_xy = data_3d[:, :, mid_z]  # XY at middle Z
    slice_xz = data_3d[:, mid_y, :]  # XZ at middle Y
    slice_yz = data_3d[mid_x, :, :]  # YZ at middle X

    if vmin is None:
        vmin = data_3d.min()
    if vmax is None:
        vmax = data_3d.max()

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # XY slice
    im1 = axes[0, 0].imshow(slice_xy.T, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(f'XY Slice (Z={mid_z})')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')
    plt.colorbar(im1, ax=axes[0, 0])

    # XZ slice
    im2 = axes[0, 1].imshow(slice_xz.T, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(f'XZ Slice (Y={mid_y})')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Z')
    plt.colorbar(im2, ax=axes[0, 1])

    # YZ slice
    im3 = axes[1, 0].imshow(slice_yz.T, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    axes[1, 0].set_title(f'YZ Slice (X={mid_x})')
    axes[1, 0].set_xlabel('Y')
    axes[1, 0].set_ylabel('Z')
    plt.colorbar(im3, ax=axes[1, 0])

    # Statistics
    axes[1, 1].axis('off')
    stats = f"""
Timestep: {timestep}
Grid: {nx} x {ny} x {nz}

Min:  {data_3d.min():.6e}
Max:  {data_3d.max():.6e}
Mean: {data_3d.mean():.6e}
Std:  {data_3d.std():.6e}
    """
    axes[1, 1].text(0.1, 0.5, stats, fontsize=10, family='monospace')

    fig.suptitle(f'Pressure Field - Timestep {timestep}')
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize BP5 pressure field data')
    parser.add_argument('nx', type=int, help='Grid X dimension')
    parser.add_argument('ny', type=int, help='Grid Y dimension')
    parser.add_argument('nz', type=int, help='Grid Z dimension')
    parser.add_argument('--file', default='snapshots.bp', help='BP5 directory path')
    parser.add_argument('--output', default='plots', help='Output directory')
    parser.add_argument('--global-scale', action='store_true', help='Use global min/max')
    parser.add_argument('--skip', type=int, default=1, help='Plot every Nth timestep')

    args = parser.parse_args()

    # Read data
    data_list, timesteps = read_bp5_file(args.file, args.nx, args.ny, args.nz)

    if data_list is None:
        print("Failed to read data")
        return 1

    # Create output dir
    os.makedirs(args.output, exist_ok=True)

    # Determine color scale
    if args.global_scale:
        vmin = min(d.min() for d in data_list)
        vmax = max(d.max() for d in data_list)
        print(f"Global scale: [{vmin:.6e}, {vmax:.6e}]")
    else:
        vmin = None
        vmax = None

    # Plot timesteps
    print(f"\nCreating plots (skip={args.skip})...")
    for i, (data, ts) in enumerate(zip(data_list, timesteps)):
        if i % args.skip != 0:
            continue

        output_file = os.path.join(args.output, f'pressure_{ts:06d}.png')
        plot_slices(data, ts, output_file, vmin=vmin, vmax=vmax)

    print(f"\nDone! Plots saved to {args.output}/")
    return 0


if __name__ == '__main__':
    sys.exit(main())
