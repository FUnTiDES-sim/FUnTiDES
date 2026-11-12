# FUnTiDES: Fast Unstructured Time Dynamic Equation Solver

FUnTiDES (Fast Unstructured Time Dynamic Equation Solver) is a high-performance computational library designed to simulate 2D and 3D wave propagation. By solving acoustic and elastic wave equations using a high-order Spectral Element Method (SEM), it handles the massive computational demands of real-world scientific simulations with exceptional accuracy. Designed primarily for use via its Python bindings, FUnTiDES combines a flexible, user-friendly Python interface with a heavily optimized, performance-portable C++ backend (powered by Kokkos) for modern High-Performance Computing (HPC) architectures.

## Included Applications

The current implementation includes proxy applications for solving the 2nd-order acoustic wave equation in 2D and 3D:

  - **SEM (Spectral Element Method)**
    A benchmark designed to simulate wave propagation using SEM, a Galerkin-based finite element method for solving partial differential equations (PDEs).

A key feature of these proxy applications is their adaptability to different programming models and HPC architectures. They are also easy to build and run, making them accessible to both researchers and developers.

## CMake Options

The following options can be used to configure your build:

| Option                     | Description                                                                        |
|----------------------------|------------------------------------------------------------------------------------|
| `COMPILE_SEM`              | Compile the Spectral Element Method solvers (default: ON)                          |
| `COMPILE_DG`               | Compile the Discontinuous Galerkin solvers (default: ON)                           |
| `COMPILE_DG_SEM`           | Compile the coupled DG-SEM solvers, forces `COMPILE_DG` and `COMPILE_SEM` (default: ON) |
| `COMPILE_DG_PADAPTIVE`     | Compile the p-adaptive DG solver, forces `COMPILE_DG` (default: ON)                |
| `COMPILE_SEM_2D`           | Compile the 2-D acoustic SEM solver (quadrilaterals in the (x, z) plane), forces `COMPILE_SEM` (default: OFF) |
| `ENABLE_PYWRAP`            | Build the Python bindings (`pyfuntides`) with pybind11 (default: OFF)              |
| `USE_MPI`                  | Compile with MPI (default: OFF)                                                    |
| `BUILD_SHARED_LIBS`        | Build shared libraries (default: ON)                                               |
| `ENABLE_COVERAGE`          | Enable code coverage. Does not work with device enabled (default: OFF)            |
| `CMAKE_INSTALL_PREFIX`     | Where to install FUnTiDES (default: `<build>/install`)                             |
| `MAX_SOLVER_ORDER`         | Max polynomial order generated for solvers (reduces compile time & binary size)    |
| `MAX_DIFFERENTIATOR_ORDER` | Max polynomial order generated for differentiators (reduces compile time)          |

The Kokkos backends (OpenMP, CUDA, ...) are not FUnTiDES options: they come from the
Kokkos installed with the TPLs (see Step 1).

**Optimizing Build Times:** Because FUnTiDES heavily relies on C++ templates for different polynomial orders (up to order 9 by default), compiling all of them can take time and memory. You can restrict the generation to lower orders if you don't need them:

```sh
cmake -S . -B build -DMAX_SOLVER_ORDER=3 -DMAX_DIFFERENTIATOR_ORDER=2
```

*Note: Advanced users can also fine-tune maximum orders per physics using specific variables like `MAX_SOLVER_ACOUSTIC_ORDER`, `MAX_SOLVER_ELASTOACOUSTIC_ORDER`, `MAX_DIFFERENTIATOR_ELASTIC_ORDER`, etc.*

## Quick Start: Build and Run

### Step 1: Install the Third-Party Libraries (TPL)

Kokkos, pybind11, the Python `kokkos` module (pykokkos-base), ADIOS2, GoogleTest and
Google Benchmark are installed by the FUnTiDES-TPL repository.

```sh
git clone --recursive https://github.com/FUnTiDES-sim/FUnTiDES-TPL.git
git clone https://github.com/FUnTiDES-sim/FUnTiDES.git
export FUNTIDES_TPL_INSTALL_DIR=$(pwd)/FUnTiDES-TPL/install
```

Then pick the target architecture.

**CPU only (OpenMP):**

```sh
FUnTiDES-TPL/install.sh --prefix=${FUNTIDES_TPL_INSTALL_DIR} --disable-cuda --enable-venv --jobs=$(nproc)
```

**NVIDIA GPU (CUDA):**

```sh
FUnTiDES-TPL/install.sh --prefix=${FUNTIDES_TPL_INSTALL_DIR} --enable-cuda --cuda-arch=89 --enable-venv --jobs=$(nproc)
```

`--cuda-arch` is the compute capability of the GPU (e.g. `80` for A100, `89` for RTX 40xx,
`90` for H100). It is printed by `nvidia-smi --query-gpu=compute_cap --format=csv`.

Notes on the TPL build:

- **Python version.** `--enable-venv` creates the virtual environment
  `${FUNTIDES_TPL_INSTALL_DIR}/tpl-venv` from the `python3` found in `PATH`. The
  `pykokkos` package requires **Python ≥ 3.11**. With an older Python, either put a newer
  one first in `PATH` before running `install.sh` (for instance `uv python install 3.12`),
  or add `--skip-pykokkos`: FUnTiDES itself only needs the `kokkos` module. Use
  `--skip-python` for a C++-only installation.
- **Memory.** Compiling pykokkos-base is memory hungry. Keep `--pykokkos-build-threads`
  low (default: 2); `--jobs` only applies to the other libraries.
- **Reruns.** Running `install.sh` again skips the components already installed;
  `--force` rebuilds everything.

For all TPL options see `FUnTiDES-TPL/install.sh --help`.

### Step 2: Set Up the Environment

`install.sh` generates a `setup_env.sh` script in the TPL install directory. Source it in
every new shell, before building or running FUnTiDES:

```sh
source ${FUNTIDES_TPL_INSTALL_DIR}/setup_env.sh
```

It sets `PATH`, `LD_LIBRARY_PATH` and `CMAKE_PREFIX_PATH`, and activates `tpl-venv` when
the TPLs were installed with `--enable-venv`.

### Step 3: Compile and Install

```sh
cd FUnTiDES
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_PYWRAP=ON \
      -DCMAKE_INSTALL_PREFIX=${FUNTIDES_TPL_INSTALL_DIR}
cmake --build build -j 8
cmake --install build
```

- FUnTiDES is installed next to the TPLs, in the same prefix: a single directory holds
  everything, and a single entry in `CMAKE_PREFIX_PATH` is enough. The files installed
  are listed in `build/install_manifest.txt` (to remove them later).
- Any other prefix works, but keep it outside the build directory, so that deleting
  `build` for a clean rebuild does not remove the installed FUnTiDES. Without
  `-DCMAKE_INSTALL_PREFIX`, FUnTiDES is installed into `build/install`.
- With CUDA, each high-order solver file can take several GB of memory to compile. If
  the build runs out of memory, use fewer jobs or lower `MAX_SOLVER_ORDER`.
- No extra option is needed for the GPU: the CUDA backend is enabled when the TPL Kokkos
  was built with CUDA.
- `-DENABLE_PYWRAP=ON` builds the Python package. CMake uses the Python of the active
  environment; it must be the one where the `kokkos` module is installed (`tpl-venv` with
  `--enable-venv`). Otherwise pass `-DPython_EXECUTABLE=/path/to/python`. Remove the
  option for a C++-only build.

### Step 4: Set Up the FUnTiDES Environment

The installation generates `funtides_env.sh` in the install directory. It sources the TPL
`setup_env.sh` (Step 2), then adds FUnTiDES to `PATH`, `LD_LIBRARY_PATH`,
`CMAKE_PREFIX_PATH` and, with `ENABLE_PYWRAP=ON`, `PYTHONPATH` (`pyfuntides` package):

```sh
source ${FUNTIDES_TPL_INSTALL_DIR}/funtides_env.sh
python -c "import kokkos, pyfuntides.model, pyfuntides.solver"
```

To load it in every new shell, add the `source` line, with the absolute path, to
`~/.bashrc`. The script can be sourced several times without duplicating entries in
FUnTiDES paths. The TPL script it relies on is found from the Kokkos used at
configuration; set `-DFUNTIDES_TPL_ENV_SCRIPT=/path/to/setup_env.sh` to choose another one,
or to an empty value to skip it.

### Using FUnTiDES from another CMake project

The installed tree provides a CMake package. Point `CMAKE_PREFIX_PATH` to the
install prefix, then:

```cmake
find_package(funtides REQUIRED)
target_link_libraries(my_app PRIVATE funtides::funtides)
```

`funtides::funtides` brings every component; a single one can be linked instead
(`funtides::model`, `funtides::solver`, `funtides::gradient`, `funtides::io`, ...).

### 2-D acoustic SEM (`COMPILE_SEM_2D`)

With `-DCOMPILE_SEM_2D=ON`, FUnTiDES also builds an acoustic SEM solver in 2-D: a mesh of
quadrilaterals in the (x, z) plane (z up, free surface on the top edge), the same
equation and time scheme as in 3-D, single process. The mesh is a `ModelUnstruct` of
dimension 2, built by `CartesianUnstruct2DBuilder` or from arrays
(`ModelUnstructData` with `dim = 2`); the solver is created with `dim = 2`:

```python
import pyfuntides.model as Model, pyfuntides.solver as Solver
params = Model.CartesianParams_f32_i32(order, ex, 1, ez, lx, 1.0, lz, False, False)  # y ignored
model = Model.CartesianUnstruct2DBuilder_f32_i32(params).get_model(True)           # free surface on top
solver = Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.UNSTRUCT,
                              Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, order, dim=2)
solver.compute_fe_init(model)
```

The data classes are the acoustic ones (`WavefieldAcoustic`, `RhsAcoustic`,
`SEMsolverDataAcoustic`) with (order+1)^2 points per element, local node
`ix + iz * (order + 1)`. Not available in 2-D yet: elastic physics, MPI, adjoint
(`update_solution_backward`), attenuation. The 3-D solvers reject a 2-D mesh. The
FUnTiDES-App repository provides the `funtides-sem2d` driver and its validation against
the exact 2-D solution; its `python/` part (makutu) has a 2-D mesh class and a 2-D
Marmousi example.

**Non-conforming blocks.** A 2-D mesh can be a set of regular blocks whose nodes are not shared
across block edges (duplicated nodes, interface nodes flagged `InteriorNode`). The blocks are glued by
a symmetric interior penalty (SIPG) applied segment by segment:
`Solver::setInterfaceOperator(face, line, dn, lag, half_inv_rho, weight, gamma, n_q)` (Python
`set_interface_operator`), one GPU thread per interface segment; the absorbing boundary skips the
interface edges. `FieldSampler` (`solver/impl/common`) sums the receivers on the device. The block
planner and the arrays given to the operator are built by makutu (`makutu.mortar`, FUnTiDES-App);
`scripts/dev/build_tree_env.sh <build dir>` sets an environment that loads a build tree before the
installed libraries (tests of a build before installing it).

### Step 5: Run Tests & Benchmarks

From the `build` directory, unit tests only:

```sh
ctest -LE "benchmark"
```

Benchmarks only, results will be stored in `build/Benchmarking` as a json file:

```sh
ctest -L benchmark
```

Or just all tests:

```sh
ctest
```

When the tests of a build directory run while another FUnTiDES is installed in a
directory listed in `LD_LIBRARY_PATH` (e.g. after sourcing `funtides_env.sh`), the test
executables load the installed libraries instead of their own (`LD_LIBRARY_PATH` comes
before their RUNPATH): run `ctest` with `LD_LIBRARY_PATH` unset in that case.

The validation against analytical solutions (`make validate`) runs the `funtides-sem`
driver and lives with it in the FUnTiDES-App repository.

### Step 6: Run Examples

The `funtides-sem` command-line driver and its YAML launcher (`run_funtides.py`) are
provided by the FUnTiDES-App repository, which builds against an installed FUnTiDES.

The Python and MPI examples are provided in the `examples` folder. They need the Python
package (Step 4) and their own requirements (`matplotlib`); the MPI example also needs
`mpi4py` (installed in `tpl-venv`):

```sh
pip install -r examples/requirements.txt

# Run a simple Python cartesian solver example
python3 examples/solver_cartesian.py --ex 100 --ey 100 --ez 100 --order 2

# Run an MPI-distributed Python cartesian solver (requires MPI)
mpirun -n 4 python3 examples/solver_cartesian_mpi.py --ex 200 --snap_interval 20
```

## 🐍 Python wrappers

### Prerequisites

The Python requirements (`numpy`, `pybind11`) are already installed in `tpl-venv`.
Without the TPL venv, install them with:

```bash
pip install -r requirements.txt
```

### Generation

Configure with `-DENABLE_PYWRAP=ON` and install (Step 3). This creates a *pyfuntides*
package in `<install prefix>/python` (`CMAKE_INSTALL_PYTHONDIR`, relative to the prefix), which contains
both the *solver* and *model* pybind modules.

### Usage

Source the generated environment script (Step 4), which makes the *pyfuntides*,
*kokkos* and *adios2* packages visible:

```bash
source <install prefix>/funtides_env.sh
```

Some examples on how to use the wrappers are available
in the [`examples`](examples/) folder.

### Tests & Benchmarks

**All commands from this section should be executed from the repository root directory\!**

To install dev python packages:

```bash
pip install -r requirements-dev.txt
```

To run python tests (default is using 6 threads):

```bash
pytest -vv -s  tests/units --kokkos-num-threads=12
pytest -vv -s  tests/integration --kokkos-num-threads=12
```

Tests marked `slow` (full simulations, such as the self-adjoint check) are skipped
unless `--run-slow` is given:

```bash
pytest -vv -s tests/integration/test_self_adjoint.py --run-slow --kokkos-num-threads=12
```

To run python benchmarks (default is using 6 threads):

```bash
pytest -vv -s tests/benchmarks --kokkos-num-threads=12
```

To run all python benchmarks (default is using 1,2,4,8,16,32,64 threads):

```bash
python scripts/benchmarks/run_pywrap_benchmarks.py --verbose
```

### Inspecting ADIOS2 Outputs

The `scripts/adios` folder provides Python helpers for BP5 outputs:

```bash
# Inspect the content of a BP5 folder
python ./scripts/adios/inspect_bp5.py <bp_folder>
python ./scripts/adios/inspect_bp5_detailed.py <bp_folder>
```

## Code Coverage

**FUnTiDES** supports code coverage analysis. To enable this feature, set the CMake option `ENABLE_COVERAGE` to `ON`.

> **Note:** The application will run significantly slower when this option is enabled. Currently, this feature is fully supported on **CPU** but has not been tested on device-specific code (GPU). Ensure that your **Kokkos** installation was not compiled with CUDA or ROCm/AMD support enabled.

### Running Tests
First, compile the code, then execute the tests using the following command:

```bash
ctest -LE "benchmark" --output-on-failure
```
This command runs all standard tests while excluding benchmarks, which can be time-consuming.

### Generating Reports
To generate a readable **HTML report**, execute the following commands:

```bash
# Capture coverage data
lcov --capture --directory . --output-file coverage.info --ignore-errors inconsistent,source

# Filter out external libraries, tests, and build files
lcov --remove coverage.info '/usr/*' '*/_deps/*' '*/tests/*' '*/buildCov/*' --output-file coverage_cleaned.info --ignore-errors inconsistent,source

# Generate the HTML report
genhtml coverage_cleaned.info \
    --output-directory html_report \
    --ignore-errors inconsistent,source,range \
    --filter missing
```
