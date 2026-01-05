# FD Solver Python Bindings

This directory contains Python bindings for the FUnTiDES Finite Difference (FD) solver, built using pybind11.

## Overview

The FD solver Python bindings provide a Python interface to the FDTD (Finite-Difference Time-Domain) acoustic wave propagation solver. This allows users to:

- Configure and run FDTD simulations from Python
- Integrate FD solver with Python-based workflows
- Prototype and test FD solver configurations interactively

## Building

The Python bindings are built as part of the main FUnTiDES build when `ENABLE_PYWRAP` is enabled:

```bash
cmake -DENABLE_PYWRAP=ON ..
make
```

This will create a Python module `fd_solver` that can be imported as:

```python
from pyfuntides import fd_solver
```

## API Overview

### Configuration Classes

#### `FdtdOptions`
Main configuration container for FDTD simulation parameters.

**Nested configuration structs:**
- `GridParams` - Grid geometry (nx, ny, nz, dx, dy, dz, mesh)
- `StencilParams` - Stencil configuration (lx, ly, lz, implem)
- `SourceParams` - Source configuration (xs, ys, zs, f0, source_order)
- `VelocityParams` - Velocity model (vmin, vmax, file_model, use_file_model)
- `TimeParams` - Time stepping (time_step, time_max, method)
- `BoundaryParams` - Boundary conditions (use_pml, use_sponge, pml_size, sponge_size, sponge_alpha)
- `OutputParams` - Output settings (save_snapshots, snapshot_interval)

**Methods:**
- `Validate()` - Validate all configuration parameters

### Component Classes

#### `FdtdGrids`
Manages computational grid geometry and model arrays.

**Methods:**
- `InitGrid(options)` - Initialize grid geometry
- `InitModelArrays(options)` - Initialize model arrays
- `nx()`, `ny()`, `nz()` - Get grid dimensions
- `dx()`, `dy()`, `dz()` - Get grid spacing

#### `FdtdStencils`
Manages finite difference stencil coefficients.

**Attributes:**
- `lx`, `ly`, `lz` - Half stencil lengths
- `coef0` - Central coefficient
- `coefx`, `coefy`, `coefz` - Stencil coefficient arrays (Kokkos Views)

**Methods:**
- `initStencilsCoefficients(options, dx, dy, dz)` - Initialize stencil coefficients

#### `FdtdKernels`
Manages computational kernels and wavefield arrays.

**Methods:**
- `initFieldsArrays(nx, ny, nz, lx, ly, lz)` - Initialize wavefield arrays

**Attributes (Kokkos Views - accessible from Python via PyKokkos):**
- `pnGlobal` - Pressure field array (2D: [n_points, 2] for time stepping)
- `RHSTerm` - Source RHS term (1D: [n_time_steps])
- `phi` - PML auxiliary field (1D: [n_points])

#### `FdtdAbcKernels`
Manages absorbing boundary conditions.

**Methods:**
- `defineSpongeBoundary(nx, ny, nz)` - Define sponge boundary conditions
- `definePML(nx, ny, nz, lx, ly, lz, dx, dy, dz, npml, dt)` - Define PML boundary conditions

**Attributes (Kokkos Views):**
- `eta` - PML damping coefficient array
- `spongeArray` - Sponge damping array

#### `FdtdSourceReceivers`
Manages source and receiver locations.

**Attributes:**
- `xsrc`, `ysrc`, `zsrc` - Source coordinates (grid indices)

### Solver Classes

#### `FdtdSolver`
Main FDTD solver class.

**Constructor:**
```python
FdtdSolver(grids, kernels, abckernels, stencils, source_receivers)
```

**Methods:**
- `compute_one_stepSB(itime, i1, i2)` - Compute one time step with sponge boundary
- `compute_one_stepPML(itime, i1, i2)` - Compute one time step with PML boundary

### Factory Functions

#### `create_fd_solver`
Factory function to create an `FdtdSolver` instance.

```python
solver = fd_solver.create_fd_solver(grids, kernels, abckernels, stencils, source_receivers)
```

#### `initialize_rhs_term`
Initialize the RHS term (source wavelet) using a Ricker wavelet.

```python
fd_solver.initialize_rhs_term(kernels, n_time_steps, dt, f0, source_order=2)
```

**Parameters:**
- `kernels` - FdtdKernels object
- `n_time_steps` - Number of time steps
- `dt` - Time step size (s)
- `f0` - Source frequency (Hz)
- `source_order` - Derivative order (0, 1, or 2), default=2

## Example Usage

```python
from pyfuntides import fd_solver

# Configure options
options = fd_solver.FdtdOptions()
options.grid.nx = 100
options.grid.ny = 100
options.grid.nz = 100
options.grid.dx = 10.0
options.grid.dy = 10.0
options.grid.dz = 10.0

options.stencil.lx = 4
options.stencil.ly = 4
options.stencil.lz = 4

options.source.xs = 50
options.source.ys = 50
options.source.zs = 50
options.source.f0 = 10.0
options.source.source_order = 2

options.boundary.use_sponge = True
options.boundary.sponge_size = 20

options.time.time_max = 0.5

# Validate configuration
options.Validate()

# Initialize components
grids = fd_solver.FdtdGrids()
grids.InitGrid(options)
grids.InitModelArrays(options)

stencils = fd_solver.FdtdStencils()
stencils.initStencilsCoefficients(options, grids.dx(), grids.dy(), grids.dz())

kernels = fd_solver.FdtdKernels()
kernels.initFieldsArrays(grids.nx(), grids.ny(), grids.nz(),
                         stencils.lx, stencils.ly, stencils.lz)

abckernels = fd_solver.FdtdAbcKernels()
abckernels.defineSpongeBoundary(grids.nx(), grids.ny(), grids.nz())

source_receivers = fd_solver.FdtdSourceReceivers()
source_receivers.xsrc = options.source.xs
source_receivers.ysrc = options.source.ys
source_receivers.zsrc = options.source.zs

# Compute time step and initialize source wavelet
vmax = options.velocity.vmax
dt = 0.5 * min(grids.dx(), grids.dy(), grids.dz()) / vmax
n_steps = int(options.time.time_max / dt)

fd_solver.initialize_rhs_term(kernels, n_steps, dt,
                               options.source.f0,
                               options.source.source_order)

# Create solver
solver = fd_solver.create_fd_solver(grids, kernels, abckernels,
                                     stencils, source_receivers)

# Run time-stepping
i1, i2 = 0, 1
for itime in range(n_steps):
    solver.compute_one_stepSB(itime, i1, i2)
    i1, i2 = i2, i1  # Swap buffers
```

For a complete working example, see:
`tests/benchmarks/python/test_fd_solver_example.py`

## Notes

- The FD solver uses Kokkos for performance portability
- Kokkos arrays (`vectorReal`, `arrayReal`) are exposed via `KokkosExp_InterOp` and can be accessed from Python using PyKokkos
- Arrays can be converted to NumPy for visualization/analysis (requires PyKokkos-base)
- Boundary condition types: Sponge or PML (Perfectly Matched Layer)
- Stencil implementations: Remez or Taylor
- Time step should satisfy CFL condition: `dt <= 0.5 * min(dx,dy,dz) / vmax`

### Working with Kokkos Arrays in Python

The FD solver exposes several Kokkos arrays that can be accessed from Python:

```python
import numpy as np
import kokkos

# After initializing kernels
kernels.initFieldsArrays(nx, ny, nz, lx, ly, lz)

# Access Kokkos arrays (read-only from bindings)
# To modify, use PyKokkos array interface
pressure = kernels.pnGlobal  # 2D array: [n_points, 2]
rhs = kernels.RHSTerm         # 1D array: [n_time_steps]

# Convert to NumPy for analysis (copy)
# pressure_np = np.array(pressure, copy=True)
```

## Related Files

- Source: `src/bindings.cpp`
- CMake: `CMakeLists.txt`
- FD Solver: `src/solver/fd/`
- Example: `tests/benchmarks/python/test_fd_solver_example.py`
