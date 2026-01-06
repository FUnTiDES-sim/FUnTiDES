#!/usr/bin/env python3
"""
Debug test for PML boundary initialization
"""

import sys
try:
    from pyfuntides import fd_solver
except ImportError as e:
    print(f"Error: Could not import fd_solver module: {e}")
    sys.exit(1)

# Initialize Kokkos
print("Initializing Kokkos...")
fd_solver.initialize_kokkos()
print(f"  Kokkos initialized: {fd_solver.is_kokkos_initialized()}")

# Create options
print("\nCreating options...")
options = fd_solver.FdtdOptions()
options.grid.nx = 50
options.grid.ny = 50
options.grid.nz = 50
options.grid.dx = 10.0
options.grid.dy = 10.0
options.grid.dz = 10.0
options.boundary.use_pml = True
options.boundary.pml_size = 10
options.stencil.lx = 2
options.Validate()
print("  Options validated")

# Initialize grid
print("\nInitializing grid...")
grids = fd_solver.FdtdGrids()
grids.InitGrid(options)
nx, ny, nz = grids.nx(), grids.ny(), grids.nz()
dx, dy, dz = grids.dx(), grids.dy(), grids.dz()
print(f"  Grid: {nx} x {ny} x {nz}")
print(f"  Spacing: {dx} x {dy} x {dz}")
print(f"  ndampx={grids.ndampx()}, ndampy={grids.ndampy()}, ndampz={grids.ndampz()}")
print(f"  x1={grids.x1()}, x2={grids.x2()}, x3={grids.x3()}")

# Initialize abckernels
print("\nInitializing ABC kernels...")
abckernels = fd_solver.FdtdAbcKernels()
print("  ABC kernels created")

# Compute dt
vmax = 4500.0
dt = 0.5 * min(dx, dy, dz) / vmax
print(f"\n  dt = {dt}, vmax = {vmax}")

# Try to initialize PML
print("\nInitializing PML...")
try:
    abckernels.definePML(grids, nx, ny, nz, dx, dy, dz, dt, vmax)
    print("  ✓ PML initialized successfully!")
except Exception as e:
    print(f"  ✗ PML initialization failed: {e}")
    import traceback
    traceback.print_exc()
    # Clean up before finalizing
    del abckernels
    del grids
    del options
    fd_solver.finalize_kokkos()
    sys.exit(1)

# Clean up Kokkos objects before finalizing
print("\nCleaning up Kokkos objects...")
del abckernels
del grids
del options

# Finalize
print("Finalizing Kokkos...")
fd_solver.finalize_kokkos()
print("  Done!")
