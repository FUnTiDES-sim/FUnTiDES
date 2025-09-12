# proxyApp: Wave simulations simplified

**proxyApp** is a collection of simplified codes that represent real scientific applications. It serves as a standard tool for evaluating and comparing the performance of various high-performance computing (HPC) systems, particularly those used for scientific simulations.

---

## Included Applications

The current implementation includes two proxy applications for solving the 2nd-order acoustic wave equation in 2D and 3D:

- **SEM (Spectral Element Method)**  
  A benchmark designed to simulate wave propagation using SEM, a Galerkin-based finite element method for solving partial differential equations (PDEs).

- **FD (Finite Difference Method)**  
  A benchmark that uses finite-difference stencil operators to simulate wave propagation and solve PDEs.

A key feature of these proxy applications is their adaptability to different programming models and HPC architectures. They are also easy to build and run, making them accessible to both researchers and developers.

---

## Supported Programming Models

The SEM proxy currently supports:

- [OpenMP](https://www.openmp.org/) — for loop-level parallelism
- [Kokkos](https://kokkos.github.io/kokkos-core-wiki/) — for performance portability

> **Note**: Kokkos is included as a Git submodule and will be compiled automatically when enabled.

---

## Supported Data Containers

The current SEM proxy supports the following data container:

- `std::vector` (default for serial and OpenMP modes)

---

## Quick Start: Build and Run

### Step 1: Compile and Install

```sh
mkdir build
cd build
cmake ..
make
```

By default, this builds the applications in sequential mode using `std::vector`. Both SEM and FD applications are compiled.

### Step 2: Run Examples

```sh
# Run SEM simulation with 100 x 100 x 100 elements
./src/main/semproxy -ex 100

# Run FD simulation
./src/main/fdproxy
```

---

## CMake Options

The following options can be used to configure your build:

| Option                 | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `COMPILE_FD`           | Enable compilation of the FD proxy (default: ON)                            |
| `COMPILE_SEM`          | Enable compilation of the SEM proxy (default: ON)                           |
| `ENABLE_CUDA`          | Enable CUDA backend (used by Kokkos)                                        |
| `ENABLE_PYWRAP`        | Enable Python bindings via pybind11 (experimental)                          |
| `USE_KOKKOS`           | Enable Kokkos support (serial by default, CUDA/OpenMP with flags)           |
| `USE_VECTOR`           | Use `std::vector` for data arrays (enabled by default unless Kokkos is used)|

---

## 🐍 Python wrappers 

### Prerequisites

The compilation of python wrappers requires _pybind11_ to be installed in your python environment.

### Generation

The proxy must be configured with `-DENABLE_PYWRAP=ON` and installed via `make install`. Optionally, you can set `-DCMAKE_INSTALL_PREFIX` to where you want to deploy the application along with the python wrappers.

This will create a _pyproxys_ package in your install directory which contains both the _solver_ and _model_ pybind modules.

```bash
(.venv) [proxys]$ ls $MY_INSTALL_DIR/pyproxys/
__init__.py  model.cpython-311-x86_64-linux-gnu.so  solver.cpython-311-x86_64-linux-gnu.so
```

This will also install _kokkos_ in your python environment, which will point to the kokkos built by the _pyproxys_ app.

```bash
(.venv) [proxys]$ ls .venv/lib/python3.11/site-packages/kokkos/
__init__.py  libpykokkos.cpython-311-x86_64-linux-gnu.so  __pycache__  pytest.ini  test  utility.py
```

### Usage

First, extend your `PYTHONPATH` to make the _pyproxys_ package visible.

```bash
export PYTHONPATH=$PYTHONPATH:$MY_INSTALL_DIR
```

Then extend your `LD_LIBRARY_PATH` so that all libraries point to the same _kokkos_ libraries that are installed in the _lib64_ folder.

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MY_INSTALL_DIR/lib64
```

There is no need to extend the `LD_LIBRARY_PATH` with the _proxys_ libraries since the python wrappers use their _RPATH_ to retrieve them in the _lib_ folder.


Some examples on how to use the wrappers are available in the [`examples`](examples/) folder.