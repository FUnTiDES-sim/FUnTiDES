#!/usr/bin/env python3
"""
Example script demonstrating the FD solver Python bindings.

This script shows how to:
1. Configure FD solver options
2. Initialize grid, stencils, kernels, and boundary conditions
3. Run a simple FDTD acoustic simulation
"""

import sys
import numpy as np
import os

try:
    from pyfuntides import fd_solver
except ImportError as e:
    print(f"Error: Could not import fd_solver module: {e}")
    print("Make sure the FD solver Python bindings are built and installed.")
    sys.exit(1)


def save_snapshot(kernels, nx, ny, nz, lx, ly, lz, itime, buffer_idx, output_dir="snapshots"):
    """Save wavefield snapshot to disk.

    Args:
        kernels: FdtdKernels object containing the wavefield
        nx, ny, nz: Grid dimensions
        lx, ly, lz: Stencil half-widths
        itime: Current time step
        buffer_idx: Current buffer index (0 or 1)
        output_dir: Directory to save snapshots
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Extract wavefield snapshot using the C++ helper function
    wavefield = fd_solver.get_wavefield_snapshot(kernels, nx, ny, nz, lx, ly, lz, buffer_idx)

    # Save as numpy binary file
    filename = os.path.join(output_dir, f"wavefield_{itime:06d}.npy")
    np.save(filename, wavefield)

    return filename


def main():
    """Run a simple FD solver example."""

    print("=" * 70)
    print("FD Solver Python Bindings Example")
    print("=" * 70)

    # =========================================================================
    # Step 0: Initialize Kokkos
    # =========================================================================
    print("\n[0] Initializing Kokkos...")
    fd_solver.initialize_kokkos()
    print(f"  ✓ Kokkos initialized (is_initialized={fd_solver.is_kokkos_initialized()})")

    # =========================================================================
    # Step 1: Configure simulation options
    # =========================================================================
    print("\n[1] Configuring simulation options...")

    options = fd_solver.FdtdOptions()

    # Grid configuration
    options.grid.nx = 200
    options.grid.ny = 200
    options.grid.nz = 200
    options.grid.dx = 10.0  # meters
    options.grid.dy = 10.0
    options.grid.dz = 10.0
    options.grid.mesh = "cartesian"

    # Stencil configuration
    options.stencil.lx = 4  # Half-width of stencil
    options.stencil.ly = 4
    options.stencil.lz = 4
    options.stencil.implem = "remez"

    # Source configuration
    options.source.xs = 50  # Source at grid center
    options.source.ys = 50
    options.source.zs = 10
    options.source.f0 = 10.0  # 10 Hz source
    options.source.source_order = 2

    # Velocity model
    options.velocity.vmin = 1500.0  # m/s
    options.velocity.vmax = 4500.0
    options.velocity.use_file_model = False

    # Time stepping
    options.time.time_step = 0.0  # Auto-compute from CFL
    options.time.time_max = 1  # 0.5 seconds
    options.time.method = "FDTD"

    # Boundary conditions
    options.boundary.use_pml = True
    options.boundary.pml_size = 2
    options.boundary.use_sponge = False
    options.boundary.sponge_size = 20
    options.boundary.sponge_alpha = 0.015

    # Output
    options.output.save_snapshots = True
    options.output.snapshot_interval = 10

    print(f"  Grid size: {options.grid.nx} x {options.grid.ny} x {options.grid.nz}")
    print(f"  Grid spacing: {options.grid.dx} x {options.grid.dy} x {options.grid.dz} m")
    print(f"  Stencil half-width: {options.stencil.lx}")
    print(f"  Source position: ({options.source.xs}, {options.source.ys}, {options.source.zs})")
    print(f"  Source frequency: {options.source.f0} Hz")
    print(f"  Boundary: Sponge (size={options.boundary.sponge_size})")

    # Validate options
    try:
        options.Validate()
        print("  ✓ Options validated successfully")
    except Exception as e:
        print(f"  ✗ Validation error: {e}")
        return 1

    # =========================================================================
    # Step 2: Initialize grid
    # =========================================================================
    print("\n[2] Initializing grid...")

    grids = fd_solver.FdtdGrids()
    grids.InitGrid(options)
    grids.InitModelArrays(options)

    nx = grids.nx()
    ny = grids.ny()
    nz = grids.nz()
    dx = grids.dx()
    dy = grids.dy()
    dz = grids.dz()

    print(f"  Grid initialized: {nx} x {ny} x {nz}")
    print(f"  Grid spacing: {dx:.2f} x {dy:.2f} x {dz:.2f} m")

    # =========================================================================
    # Step 3: Initialize stencils
    # =========================================================================
    print("\n[3] Initializing stencils...")

    stencils = fd_solver.FdtdStencils()
    stencils.initStencilsCoefficients(options, dx, dy, dz)

    print(f"  Stencil half-widths: lx={stencils.lx}, ly={stencils.ly}, lz={stencils.lz}")
    print(f"  Central coefficient: {stencils.coef0:.6f}")

    # =========================================================================
    # Step 4: Initialize kernels
    # =========================================================================
    print("\n[4] Initializing kernels...")

    kernels = fd_solver.FdtdKernels()
    kernels.initFieldsArrays(nx, ny, nz, stencils.lx, stencils.ly, stencils.lz)

    print("  ✓ Wavefield arrays allocated")

    # =========================================================================
    # Step 5: Initialize boundary conditions
    # =========================================================================
    print("\n[5] Initializing boundary conditions...")

    abckernels = fd_solver.FdtdAbcKernels()

    if options.boundary.use_sponge:
        abckernels.defineSpongeBoundary(nx, ny, nz)
        print(f"  ✓ Sponge boundary initialized (size={options.boundary.sponge_size})")
    elif options.boundary.use_pml:
        # Compute time step (simple CFL)
        vmax = options.velocity.vmax
        dt = 0.5 * min(dx, dy, dz) / vmax
        abckernels.definePML(grids, nx, ny, nz, dx, dy, dz, dt, vmax)
        print(f"  ✓ PML boundary initialized (size={options.boundary.pml_size})")

    # =========================================================================
    # Step 6: Initialize source/receivers
    # =========================================================================
    print("\n[6] Initializing source...")

    source_receivers = fd_solver.FdtdSourceReceivers()
    source_receivers.xsrc = options.source.xs
    source_receivers.ysrc = options.source.ys
    source_receivers.zsrc = options.source.zs

    print(f"  Source location: ({source_receivers.xsrc}, {source_receivers.ysrc}, {source_receivers.zsrc})")

    # =========================================================================
    # Step 7: Initialize source wavelet (RHS term)
    # =========================================================================
    print("\n[7] Initializing source wavelet...")

    # Compute time step using CFL condition
    vmax = options.velocity.vmax
    dt = 0.5 * min(dx, dy, dz) / vmax
    n_steps = int(options.time.time_max / dt)

    print(f"  Time step: {dt:.6f} s (CFL limited)")
    print(f"  Number of steps: {n_steps}")

    # Initialize Ricker wavelet
    fd_solver.initialize_rhs_term(kernels, n_steps, dt,
                                   options.source.f0,
                                   options.source.source_order)

    print(f"  ✓ Ricker wavelet initialized (f0={options.source.f0} Hz, order={options.source.source_order})")

    # =========================================================================
    # Step 8: Create solver
    # =========================================================================
    print("\n[8] Creating FD solver...")

    solver = fd_solver.create_fd_solver(grids, kernels, abckernels, stencils, source_receivers)

    print("  ✓ FD solver created successfully")

    # =========================================================================
    # Step 9: Run time-stepping loop (demonstration)
    # =========================================================================
    print("\n[9] Running time-stepping demonstration...")

    print(f"  Total simulation time: {n_steps * dt:.6f} s")
    print(f"  Using {'sponge' if options.boundary.use_sponge else 'PML'} boundary conditions")

    # Run a few time steps as demonstration
    n_demo_steps = min(1000, n_steps)
    print(f"\n  Running {n_demo_steps} time steps as demonstration...")

    # Initialize snapshot counter
    n_snapshots_saved = 0

    i1, i2 = 0, 1
    for itime in range(n_demo_steps):
        if options.boundary.use_sponge:
            solver.compute_one_stepSB(itime, i1, i2)
        else:
            solver.compute_one_stepPML(itime, i1, i2)
        i1, i2 = i2, i1  # Swap buffers

        # Save snapshot if enabled and at the right interval
        if options.output.save_snapshots and (itime % options.output.snapshot_interval == 0):
            filename = save_snapshot(kernels, nx, ny, nz, stencils.lx, stencils.ly, stencils.lz, itime, i2)
            n_snapshots_saved += 1
            if n_snapshots_saved <= 3:  # Only print first few to avoid spam
                print(f"    Saved snapshot: {filename}")

        if (itime + 1) % 5 == 0:
            print(f"    Completed time step {itime + 1}/{n_demo_steps}")

    print(f"  ✓ Time-stepping completed successfully")

    if options.output.save_snapshots:
        print(f"  ✓ Saved {n_snapshots_saved} snapshots to 'snapshots/' directory")

    # Note: For a full simulation, continue the loop and add output routines
    print("\n  Note: Full simulation would continue for all", n_steps, "time steps")
    print("  and include output routines for visualization.")

    print("\n" + "=" * 70)
    print("FD Solver Python bindings example completed successfully!")
    print("=" * 70)

    # =========================================================================
    # Step 10: Clean up and finalize Kokkos
    # =========================================================================
    print("\n[10] Cleaning up Kokkos objects...")
    del solver
    del source_receivers
    del abckernels
    del kernels
    del stencils
    del grids
    del options

    print("\n[11] Finalizing Kokkos...")
    fd_solver.finalize_kokkos()
    print("  ✓ Kokkos finalized")

    return 0


if __name__ == "__main__":
    sys.exit(main())
