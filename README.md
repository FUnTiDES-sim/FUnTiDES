# FUnTiDES: Fast Unstructured Time Dynamic Equation Solver

FUnTiDES (Fast Unstructured Time Dynamic Equation Solver) is a high-performance computational library designed to simulate 2D and 3D wave propagation. By solving acoustic and elastic wave equations using a high-order Spectral Element Method (SEM), it handles the massive computational demands of real-world scientific simulations with exceptional accuracy. Designed primarily for use via its Python bindings, FUnTiDES combines a flexible, user-friendly Python interface with a heavily optimized, performance-portable C++ backend (powered by Kokkos) for modern High-Performance Computing (HPC) architectures.

---

## Included Applications

The current implementation includes two proxy applications for solving the 2nd-order acoustic wave equation in 2D and 3D:

- **SEM (Spectral Element Method)**
  A benchmark designed to simulate wave propagation using SEM, a Galerkin-based finite element method for solving partial differential equations (PDEs).

- **FD (Finite Difference Method)**
  A benchmark that uses finite-difference stencil operators to simulate wave propagation and solve PDEs.

A key feature of these proxy applications is their adaptability to different programming models and HPC architectures. They are also easy to build and run, making them accessible to both researchers and developers.

---


## Quick Start: Build and Run

### Step 1: Install Third-Party Libraries (TPL)

Before compiling the main project, you must install the required dependencies via the FUnTiDES-TPL repository. 

```sh
# Clone the TPL repository
git clone --recurse-submodules https://github.com/FUnTiDES-sim/FUnTiDES-TPL.git

# Install TPL for CPU
export FUNTIDES_TPL_INSTALL_DIR=$(pwd)/FUnTiDES-TPL/install
FUnTiDES-TPL/install.sh --prefix=${FUNTIDES_TPL_INSTALL_DIR} --disable-cuda --disable-mpi --no-venv --jobs=$(nproc)
````

For more TPL builds options see `FUnTiDES-TPL/install.sh --help`.

### Step 2: Compile and Install

Once the TPLs are installed and your environment is configured, you can build the applications.

```sh
source FUnTiDES-TPL/setup_tpl_env.sh ${FUNTIDES_TPL_INSTALL_DIR}
cd FUnTiDES
mkdir build
cd build
cmake ..
make install
```

### Step 3: Run Tests & Benchmarks

Unit tests only:

```sh
ctest -LE "benchmark|validation"
```

Benchmarks only, results will be stored in `build/Benchmarking` as a json file:

```sh
ctest -L benchmark
```

Validation tests (compare against analytical solutions):

```sh
# All 24 tests (2 mesh types × 3 orders × 2 physics × 2 model-on-nodes)
make validate

# Quick validation (12 tests, without model-on-nodes)
make validate_no_model_on_nodes

# Specific groups
make validate_acoustic      # Acoustic tests only
make validate_cartesian     # Cartesian mesh tests only

# Fine-grained filtering
ctest -L validation_acoustic        # Acoustic tests
ctest -L order2                     # Order 2 tests
ctest -L mesh_ucartesian            # Ucartesian mesh tests
ctest -R acoustic_order2_cartesian  # Specific test
```

Or just all tests:

```sh
ctest
```

> **Note**: Validation tests require analytical reference solutions in `tests/analyticalsolution/` (P.dat for acoustic, Ux.dat for elastic). You also need to install the adios2 python module to run the validation tests as it runs a python script for receiver output.

### Step 4: Run Examples

Inside the `build/bin` folder, you can run the generated executable binaries. For the SEM proxy, use `semproxy` (or `funtides-sem`) and specify your mesh and implementation preferences.

```sh
# Run SEM simulation (e.g., Acoustic, Cartesian Mesh, Makutu Implementation, Order 2)
./funtides-sem --ex 100 --ey 100 --ez 100 --method=sem --implem=makutu --mesh=cartesian -o 2 --dt 0.001 --timemax 1.5

# Run FD simulation
./funtidesfd

# Run validation with custom parameters
./validate_solution --order 2 --mesh ucartesian --elastic --is-model-on-nodes
```

You can also run the Python and MPI examples provided in the `examples/fe` folder:

```sh
# Run a simple Python cartesian solver example
python3 examples/fe/solver_cartesian.py --ex 100 --ey 100 --ez 100 --order 2

# Run an MPI-distributed Python cartesian solver (requires MPI)
mpirun -n 4 python3 examples/fe/solver_cartesian_mpi.py --ex 200 --snap_interval 20
```

> **Note**: A dedicated python env can be made via TPLs.

-----

## CMake Options

The following options can be used to configure your build:

| Option          | Description                                        |
|-----------------|----------------------------------------------------|
| `COMPILE_FD`    | Enable compilation of the FD proxy (default: ON)   |
| `COMPILE_SEM`   | Enable compilation of the SEM proxy (default: ON)  |
| `ENABLE_PYWRAP` | Enable Python bindings via pybind11 (experimental) |

-----

## 🐍 Python wrappers

### Prerequisites

To install python requirements:

```bash
pip install -r requirements.txt
```

### Generation

The proxy must be configured with `-DENABLE_PYWRAP=ON` and installed via `make install`. Optionally, you can set `-DCMAKE_INSTALL_PREFIX` to where you want to deploy the application along with the python wrappers.

This will create a *pyfuntides* package in your install directory which contains both the *solver* and *model* pybind modules.

```bash
(.venv) [proxys]$ ls $MY_INSTALL_DIR/python/pyfuntides/
__init__.py  model.cpython-311-x86_64-linux-gnu.so  solver.cpython-311-x86_64-linux-gnu.so
```

This will also install *kokkos* in your python environment, which will point to the kokkos built by the *pyfuntides* app.

```bash
(.venv) [proxys]$ ls .venv/lib/python3.11/site-packages/kokkos/
__init__.py  libpykokkos.cpython-311-x86_64-linux-gnu.so  __pycache__  pytest.ini  test  utility.py
```

If you do not have write access on your python environment, it will install it under *$MY\_INSTALL\_DIR/lib/python3.11/site-packages/kokkos*. In that case you will have to extend your python path with this directory.

### Usage

First, extend your `PYTHONPATH` to make the *pyfuntides* and *adios* packages visible.

```bash
export PYTHONPATH=$PYTHONPATH:$MY_INSTALL_DIR/python
```

If needed (kokkos could not write in your python environment), also extend your `PYTHONPATH` to make the kokkos package visible.

```bash
export PYTHONPATH=$PYTHONPATH:$MY_INSTALL_DIR/lib/python3.11/site-packages
```

Then extend your `LD_LIBRARY_PATH` so that all libraries point to the same *kokkos* libraries that are installed in the *lib64* folder.

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MY_INSTALL_DIR/lib64
```

There is no need to extend the `LD_LIBRARY_PATH` with the *proxys* libraries since the python wrappers use their *RPATH* to retrieve them in the *lib* folder.

Some examples on how to use the wrappers are available in the [`examples`](examples/) folder.

### Tests & Benchmarks

**All commands from this section should be executed from the repository root directory\!**

To install dev python packages:

```bash
pip install -r requirements-dev.txt
```

To run basic python unit tests (default is using 6 threads):

```bash
pytest -vv -s  tests/units --kokkos-num-threads=12
```

To run basic python integration tests (default is using 6 threads)
```bash
pytest -vv -s  tests/integration --kokkos-num-threads=12
```

To run python benchmarks (default is using 6 threads):

```bash
pytest -vv -s tests/benchmarks --kokkos-num-threads=12
```

To run all python benchmarks (default is using 1,2,4,8,16,32,64 threads):

```bash
python scripts/benchmarks/run_pywrap_benchmarks.py --verbose
```

### Plotting Receivers and Snapshots

To plot the snapshots we provide a python script:

```bash
python ./scripts/adios/adios_cartesian_snap_viz.py 201 201 201 --file snapshots.bp --slice
```

where 201 values should be replaced by number of nodes on x, y, and z. And file corresponds to the `snapshots.bp` folder with bp5 files.

For the receivers:

```bash
python ./scripts/adios/adios_single_receiver_viz.py
```

within the folder containing the `receivers.bp` folder.
