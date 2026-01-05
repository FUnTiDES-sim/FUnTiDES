#!/usr/bin/env python3
"""Simple test to verify FD solver Python bindings work."""

import sys

try:
    from pyfuntides import fd_solver
except ImportError as e:
    print(f"Error: Could not import fd_solver module: {e}")
    sys.exit(1)

print("Testing FD Solver Python Bindings")
print("=" * 50)

# Initialize Kokkos
print("\n1. Initializing Kokkos...")
fd_solver.initialize_kokkos()
print(f"   Kokkos initialized: {fd_solver.is_kokkos_initialized()}")

# Test configuration
print("\n2. Creating FdtdOptions...")
options = fd_solver.FdtdOptions()
options.grid.nx = 50
options.grid.ny = 50
options.grid.nz = 50
options.grid.dx = 10.0
options.grid.dy = 10.0
options.grid.dz = 10.0
options.stencil.lx = 2
options.stencil.ly = 2
options.stencil.lz = 2
options.source.xs = 25
options.source.ys = 25
options.source.zs = 25
options.source.f0 = 10.0
options.velocity.vmin = 1500.0
options.velocity.vmax = 4500.0
options.boundary.use_pml = False
options.boundary.use_sponge = True
options.time.time_max = 0.1
print(f"   Grid: {options.grid.nx} x {options.grid.ny} x {options.grid.nz}")

# Validate
print("\n3. Validating options...")
try:
    options.Validate()
    print("   ✓ Validation passed")
except Exception as e:
    print(f"   ✗ Validation failed: {e}")
    sys.exit(1)

# Initialize grid
print("\n4. Initializing grid...")
grids = fd_solver.FdtdGrids()
grids.InitGrid(options)
grids.InitModelArrays(options)
print(f"   ✓ Grid: {grids.nx()} x {grids.ny()} x {grids.nz()}")

# Initialize stencils
print("\n5. Initializing stencils...")
stencils = fd_solver.FdtdStencils()
stencils.initStencilsCoefficients(options, grids.dx(), grids.dy(), grids.dz())
print(f"   ✓ Stencil lengths: {stencils.lx}, {stencils.ly}, {stencils.lz}")

# Initialize kernels
print("\n6. Initializing kernels...")
kernels = fd_solver.FdtdKernels()
kernels.initFieldsArrays(grids.nx(), grids.ny(), grids.nz(),
                         stencils.lx, stencils.ly, stencils.lz)
print("   ✓ Kernels initialized")

# Initialize source wavelet
print("\n7. Initializing source wavelet...")
dt = 0.5 * min(grids.dx(), grids.dy(), grids.dz()) / options.velocity.vmax
n_steps = int(options.time.time_max / dt)
fd_solver.initialize_rhs_term(kernels, n_steps, dt, options.source.f0, 2)
print(f"   ✓ Source initialized ({n_steps} time steps, dt={dt:.6f}s)")

print("\n" + "=" * 50)
print("✓ All basic tests passed!")
print("=" * 50)

# Finalize
fd_solver.finalize_kokkos()
print("\nKokkos finalized")
