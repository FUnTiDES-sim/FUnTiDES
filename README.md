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
#### Using JSON Configuration (Recommended)

The SEM proxy now supports JSON configuration files for easier parameter management:
```sh
# Create a configuration file
cat > config.json << 'EOF'
{
  "simulation": {
    "order": 4,
    "method": "sem",
    "implementation": "makutu",
    "mesh": "cartesian",
    "dt": 0.006,
    "timemax": 1.5,
    "autodt": false
  },
  "domain": {
    "ex": 100,
    "ey": 100,
    "ez": 100,
    "lx": 2000.0,
    "ly": 2000.0,
    "lz": 2000.0
  },
  "source": {
    "x": 1010.0,
    "y": 1010.0,
    "z": 1010.0
  },
  "receiver": {
    "x": 1410.0,
    "y": 1010.0,
    "z": 1010.0
  },
  "snapshots": {
    "enabled": true,
    "time_interval": 10
  },
  "boundaries": {
    "size": 200.0,
    "surface_sponge": true,
    "taper_delta": 0.015
  },
  "model": {
    "is_on_nodes": false,
    "is_elastic": true
  }
}
EOF

# Run with JSON configuration
./bin/semproxy -i config.json
```

#### JSON Configuration Reference

All fields are optional; omitted values use defaults:

| JSON Field | Description | Default |
|------------|-------------|---------|
| `simulation.order` | Polynomial order | 2 |
| `simulation.method` | Method: `sem` or `dg` | `sem` |
| `simulation.implementation` | Implementation: `makutu` or `shiva` | `makutu` |
| `simulation.mesh` | Mesh type: `cartesian` or `ucartesian` | `cartesian` |
| `simulation.dt` | Time step (seconds) | 0.006 |
| `simulation.timemax` | Simulation duration (seconds) | 0.7 |
| `simulation.autodt` | Auto-compute dt via CFL | false |
| `domain.ex`, `domain.ey`, `domain.ez` | Number of elements | 50, 50, 50 |
| `domain.lx`, `domain.ly`, `domain.lz` | Domain size (meters) | 2000, 2000, 2000 |
| `source.x`, `source.y`, `source.z` | Source position (meters) | 1010, 1010, 1010 |
| `receiver.x`, `receiver.y`, `receiver.z` | Receiver position (meters) | 1410, 1010, 1010 |
| `snapshots.enabled` | Enable snapshots | false |
| `snapshots.time_interval` | Snapshot interval (iterations) | 10 |
| `boundaries.size` | Sponge boundary size (meters) | 0.0 |
| `boundaries.surface_sponge` | Surface as non-sponge | false |
| `boundaries.taper_delta` | Taper delta coefficient | 0.015 |
| `model.is_on_nodes` | Model defined on nodes vs elements | false |
| `model.is_elastic` | Elastic vs acoustic simulation | false |

#### Minimal JSON Example
```json
{
  "simulation": {
    "order": 6
  },
  "domain": {
    "ex": 100,
    "ey": 100,
    "ez": 100
  }
}
```

All other parameters will use default values.

#### Legacy CLI Mode (Still Supported)

Command-line options are still available for backward compatibility:
```sh
# Run SEM simulation with CLI options (legacy mode)
./bin/semproxy --order 4 --ex 100 --ey 100 --ez 100 --lx 2000 --ly 2000 --lz 2000

# Run FD simulation
./bin/fdproxy
```

For a full list of CLI options:
```sh
./bin/semproxy --help
```

> **Note**: When using `-i <config.json>`, CLI options are ignored. Use either JSON or CLI, not both.


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
