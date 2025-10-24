# FUnTiDES: Fast Unstructured Time Dynamic Equation Solver

**FUnTiDES** is a collection of simplified codes that represent real scientific applications. It serves as a standard tool for evaluating and comparing the performance of various high-performance computing (HPC) systems, particularly those used for scientific simulations.

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

- [Kokkos](https://kokkos.github.io/kokkos-core-wiki/) — for performance portability

> **Note**: Kokkos is included as a Git submodule and will be compiled automatically when enabled.

---

## Supported Data Containers

The current SEM proxy supports the following data container:

- `std::vector` (default for serial )

---

## Quick Start: Build and Run

### Step 1: Compile and Install

```sh
mkdir build
cd build
cmake ..
make install
```

By default, this builds the applications in sequential mode using `std::vector`.
Both SEM and FD applications are compiled.

### Step 2: Run Tests & Benchmarks

Unit tests only
```sh
ctest -LE benchmark
```

Benchmarks only, results will be stored in results generated in `build/Benchmarking` as a json file.
```sh
ctest -L benchmark
```

Or just both
```sh
ctest
```

### Step 3: Run Examples

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

To install python requirements
```bash
pip install -r requirements.txt
```

### Generation

The proxy must be configured with `-DENABLE_PYWRAP=ON` and installed via `make install`. Optionally, you can set `-DCMAKE_INSTALL_PREFIX` to where you want to deploy the application along with the python wrappers.

This will create a _pyfuntides_ package in your install directory which contains both the _solver_ and _model_ pybind modules.

```bash
(.venv) [proxys]$ ls $MY_INSTALL_DIR/python/pyfuntides/
__init__.py  model.cpython-311-x86_64-linux-gnu.so  solver.cpython-311-x86_64-linux-gnu.so
```

This will also install _kokkos_ in your python environment, which will point to the kokkos built by the _pyfuntides_ app.

```bash
(.venv) [proxys]$ ls .venv/lib/python3.11/site-packages/kokkos/
__init__.py  libpykokkos.cpython-311-x86_64-linux-gnu.so  __pycache__  pytest.ini  test  utility.py
```

If you do not have write access on your python environment, it will install it under _$MY_INSTALL_DIR/lib/python3.11/site-packages/kokkos_.
In that case you will have to extend your python path with this directory.

### Usage

First, extend your `PYTHONPATH` to make the _pyfuntides_ and _adios_ package visible.

```bash
export PYTHONPATH=$PYTHONPATH:$MY_INSTALL_DIR/python
```

If needed (kokkos could not write in your python environment), also extend your `PYTHONPATH` to make the kokkos package visible.
```bash
export PYTHONPATH=$PYTHONPATH:$MY_INSTALL_DIR/lib/python3.11/site-packages
```

Then extend your `LD_LIBRARY_PATH` so that all libraries point to the same _kokkos_ libraries that are installed in the _lib64_ folder.

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MY_INSTALL_DIR/lib64
```

There is no need to extend the `LD_LIBRARY_PATH` with the _proxys_ libraries since the python wrappers use their _RPATH_ to retrieve them in the _lib_ folder.


Some examples on how to use the wrappers are available in the [`examples`](examples/) folder.

### Tests & Benchmarks

**All commands from this section should be executed from the repository root directory!**

To install dev python packages
```bash
pip install -r requirements-dev.txt
```

To run basic python unit tests (default is using 6 threads)
```bash
pytest -vv -s  tests/units
```

To run basic python unit tests with more threads
```bash
pytest -vv -s  tests/units --threads 12
```

To run python benchmarks (default is using 6 threads)
```bash
pytest -vv -s tests/benchmarks/python
```

To run python benchmarks with more threads
```bash
pytest -vv -s tests/benchmarks/python --threads 12
```

To run all python benchmarks (default is using 1,2,4,8,16,32,64 threads)
```bash
python scripts/benchmarks/run_pywrap_benchmarks.py --verbose
```

### Ploting Receivers and Snapshots

To plot the snapshots we provide a python script:
```bash
python ./scripts/adios/adios_cartesian_snap_viz.py 201 201 201 --file snapshots.bp --slice
```
where 201 values should be replaced by number of nodes on x y and z. And file correpond to the `snapshots.bp` folder with bp5 files.

For the receivers:
``` bash
python ./scripts/adios/adios_single_receiver_viz.py
```
within the folder containing the `receivers.bp` folder.
