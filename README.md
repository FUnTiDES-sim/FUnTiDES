# 🌊 FUnTiDES

> **F**ast **Un**structured **Ti**me **D**ynamic **E**quation **S**olver

A high-performance framework for seismic and acoustic wave propagation simulations, designed for modern HPC architectures.

[![C++17](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.12+-064F8C.svg)](https://cmake.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## 🎯 Overview

FUnTiDES provides production-ready proxy applications for solving the 2nd-order acoustic wave equation in 2D and 3D. It serves as a standard benchmark for evaluating HPC system performance and supports multiple numerical methods and parallel programming models.

### Key Features

- 🚀 **Multiple Numerical Methods**: Spectral Element Method (SEM) and Finite Difference Method (FD)
- 🔧 **Flexible Parallelism**: OpenMP, Kokkos (CUDA, HIP, SYCL)
- 🐍 **Python Bindings**: Experimental PyBind11 wrappers
- 📊 **Advanced I/O**: ADIOS2 support for snapshot and receiver output
- ⚡ **Performance Portable**: Single codebase, multiple architectures

---

## 🚀 Quick Start

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.12 or higher
- (Optional) CUDA Toolkit for GPU support
- (Optional) Python 3.7+ for Python bindings

### Build

```bash
# Clone the repository
git clone --recursive https://github.com/your-org/funtides.git
cd funtides

# Configure and build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Optional: Install
make install
```

### Run Your First Simulation

```bash
# Spectral Element Method - 100³ element mesh
./src/main/semproxy -ex 100

# Finite Difference Method - default configuration
./src/main/fdproxy
```

---

## 📦 Numerical Methods

### Spectral Element Method (SEM)

High-order Galerkin finite element method ideal for complex geometries and accurate wave propagation in heterogeneous media.

**Use when:**
- High accuracy is required
- Complex domain geometries
- Heterogeneous material properties

### Finite Difference Method (FD)

Efficient stencil-based approach for regular grids, optimized for large-scale 3D simulations.

**Use when:**
- Large computational domains
- Regular Cartesian meshes
- Memory bandwidth is critical

---

## ⚙️ Configuration Options

### CMake Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `COMPILE_FD` | `ON` | Build Finite Difference proxy |
| `COMPILE_SEM` | `ON` | Build Spectral Element proxy |
| `USE_KOKKOS` | `OFF` | Enable Kokkos for performance portability |
| `ENABLE_CUDA` | `OFF` | Enable NVIDIA GPU support via Kokkos |
| `ENABLE_PYWRAP` | `OFF` | Build Python bindings (experimental) |
| `USE_VECTOR` | `ON` | Use `std::vector` containers |

### Build Examples

**OpenMP parallel:**
```bash
cmake -DUSE_OPENMP=ON ..
```

**CUDA GPU acceleration:**
```bash
cmake -DUSE_KOKKOS=ON -DENABLE_CUDA=ON ..
```

**Python bindings:**
```bash
cmake -DENABLE_PYWRAP=ON -DCMAKE_INSTALL_PREFIX=/path/to/install ..
```

---

## 🎮 Usage Guide

### SEM Proxy

#### Command Line Interface

```bash
semproxy [OPTIONS]
```

#### Core Options

##### Mesh Configuration
```bash
--ex 100 --ey 100 --ez 100     # Elements per direction (Cartesian)
--lx 1000 --ly 1000 --lz 1000  # Domain size in meters
--mesh cartesian               # Mesh type: cartesian|ucartesian
```

##### Discretization
```bash
-o 4                           # Polynomial order (spectral accuracy)
--method sem                   # Method: sem|dg (Discontinuous Galerkin)
--implem makutu                # Implementation: makutu|shiva
```

##### Time Stepping
```bash
--dt 0.001                     # Time step (seconds)
--timemax 2.0                  # Simulation duration (seconds)
--auto-dt                      # Auto-compute from CFL condition
```

##### Boundary Conditions
```bash
--boundaries-size 100          # Absorbing boundary thickness (meters)
--sponge-surface               # Exclude surface nodes from sponge
--taper-delta 0.95             # Sponge layer damping coefficient
```

##### Output
```bash
-s                             # Enable snapshots
--snap-interval 50             # Snapshots every N iterations
```

##### Model Configuration
```bash
--is-model-on-nodes            # Model defined on nodes vs elements
```

#### Example Workflows

```bash
# Standard 3D simulation with auto time-stepping
semproxy -o 4 --ex 50 --ey 50 --ez 50 \
         --lx 5000 --ly 5000 --lz 5000 \
         --auto-dt --timemax 3.0 \
         -s --snap-interval 100

# High-order 2D simulation with sponge boundaries
semproxy -o 6 --ex 200 --ey 200 --ez 1 \
         --method sem --implem makutu \
         --boundaries-size 200 --sponge-surface \
         --dt 0.0005 --timemax 1.5
```

---

### FD Proxy

#### Command Line Interface

```bash
fdproxy [OPTIONS]
```

#### Core Options

##### Grid Configuration
```bash
--nx 201 --ny 201 --nz 201    # Grid dimensions (nodes)
--dx 10 --dy 10 --dz 10        # Grid spacing (meters)
--mesh cartesian               # Mesh type
```

##### Stencil Configuration
```bash
--lx 4 --ly 4 --lz 4           # Stencil half-width per direction
--implem remez                 # Implementation: remez|taylor
```

##### Source Configuration
```bash
--f0 25.0                      # Peak frequency (Hz)
--xs 1000 --ys 1000 --zs 1000  # Source position (meters, -1=center)
--sourceOrder 2                # Time derivative order
```

##### Velocity Model
```bash
--vmin 1500 --vmax 3500        # Velocity range (m/s)
--fileModel model.bin          # Load velocity model from file
```

##### Time Stepping
```bash
--timeStep 0.001               # Time step (seconds, 0=auto CFL)
--timeMax 2.0                  # Simulation duration (seconds)
--method FDTD                  # Time stepping method
```

##### Boundary Conditions
```bash
--usePML                       # Enable PML absorbing boundaries
--pmlSize 20                   # PML thickness (grid points)
--spongeSize 10                # Sponge layer thickness
--spongeAlpha 0.9              # Sponge damping coefficient
```

##### Output
```bash
--saveSnapShots                # Save wavefield snapshots
--snapShotInterval 100         # Snapshot frequency (time steps)
```

#### Example Workflows

```bash
# Small acoustic simulation with PML boundaries
fdproxy --nx 101 --ny 101 --nz 101 \
        --dx 10 --dy 10 --dz 10 \
        --f0 20 --timeMax 1.0 \
        --usePML --pmlSize 15 \
        --saveSnapShots --snapShotInterval 50

# High-frequency simulation with custom velocity
fdproxy --nx 301 --ny 301 --nz 301 \
        --vmin 2000 --vmax 4500 \
        --f0 50 --sourceOrder 2 \
        --timeStep 0 --timeMax 3.0

# 2D simulation with Remez stencil
fdproxy --nx 501 --ny 501 --nz 1 \
        --dx 5 --dy 5 \
        --lx 6 --ly 6 --implem remez \
        --f0 30 --xs -1 --ys -1 \
        --timeStep 0 --timeMax 2.5
```

---

## 🐍 Python Interface

### Installation

```bash
# Build with Python support
cmake -DENABLE_PYWRAP=ON -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
make install

# Configure environment
export PYTHONPATH=$HOME/.local:$PYTHONPATH
export LD_LIBRARY_PATH=$HOME/.local/lib64:$LD_LIBRARY_PATH
```

### Quick Example

```python
from pyproxys import solver, model

# Create simulation
sim = solver.FDSimulation(nx=201, ny=201, nz=201)
sim.set_source(x=1000, y=1000, z=1000, f0=25.0)

# Run simulation
sim.run(time_max=2.0)

# Access results
snapshots = sim.get_snapshots()
receivers = sim.get_receivers()
```

### Testing

```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run unit tests
pytest -vv tests/units

# Run benchmarks
pytest -vv tests/benchmarks

# Generate benchmark plots
pytest --benchmark-histogram=plot tests/benchmarks
```

---

## 📊 Visualization

### Snapshot Visualization

```bash
# Visualize 3D snapshots (with slicing)
python scripts/adios/adios_cartesian_snap_viz.py 201 201 201 \
       --file snapshots.bp --slice
```

### Receiver Traces

```bash
# Plot receiver time series
cd output_directory
python scripts/adios/adios_single_receiver_viz.py
```

---

## 🏗️ Project Structure

```
funtides/
├── src/
│   ├── main/           # Proxy application entry points
│   ├── sem/            # Spectral Element Method implementation
│   ├── fd/             # Finite Difference implementation
│   └── common/         # Shared utilities and kernels
├── examples/           # Python usage examples
├── scripts/            # Visualization and analysis tools
├── tests/
│   ├── units/          # Unit tests
│   └── benchmarks/     # Performance benchmarks
└── external/           # Third-party dependencies (Kokkos)
```

---

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Code Style

- C++: Google C++ Style Guide
- Python: PEP 8
- Editor: Doom Emacs configuration included

---

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgments

- Kokkos team for the performance portability framework
- ADIOS2 for parallel I/O capabilities
- The scientific computing community

---

## 📚 References

If you use FUnTiDES in your research, please cite:

```bibtex
@software{funtides2024,
  title = {FUnTiDES: Fast Unstructured Time Dynamic Equation Solver},
  author = {Your Team},
  year = {2024},
  url = {https://github.com/your-org/funtides}
}
```

---
