#!/bin/bash
# Build the solver_tests suite (non-O3, sm_103) and run the solver unit tests
# on elba-b02. Part of the FUnTiDES src/solver refactor validation.
#
# Usage: run_solver_tests.sh [test-name-filter]
#   default filter: all tests under tests/units/solver
set -euo pipefail

NODE=${SEMPROXY_BUILD_NODE:-elba-b02}
PARTITION=${SEMPROXY_BUILD_PARTITION:-COMPUTE_GPU_B300}
if [ "$(hostname -s)" != "$NODE" ]; then
  exec /usr/bin/srun --partition="$PARTITION" --nodelist="$NODE" \
    --cpus-per-task=32 --mem=64G --time=04:00:00 --quiet \
    --comment="FUnTiDES,LLM4CSE,OT/R&D/SUB/SC" \
    /appli_RD/j1039830/FUnTiDES/scripts/run_solver_tests.sh "$@"
fi

export PATH=/appli_RD/j1039830/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/appli_RD/j1039830/local/cuda/lib64:/appli_RD/j1039830/local/lib64:${LD_LIBRARY_PATH:-}

CMAKE=/appli_RD/j1039830/local/cmake/bin/cmake
CTEST=/appli_RD/j1039830/local/cmake/bin/ctest
SRC=/appli_RD/j1039830/FUnTiDES
BUILD=/appli_RD/j1039830/FUnTiDES/build

# Configure with the test suite enabled (idempotent).
"${CMAKE}" -S "${SRC}" -B "${BUILD}" \
  -DCMAKE_BUILD_TYPE= \
  -DPython_EXECUTABLE=/appli_RD/j1039830/local/tpl-venv/bin/python \
  -DKokkos_DIR=/appli_RD/j1039830/local/lib64/cmake/Kokkos

# Build only the solver test targets.
"${CMAKE}" --build "${BUILD}" --target solver_tests -j32

# Pick a free GPU for the test run (gtest_discover_tests OOMs on busy GPUs).
FREE_GPU=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
  | sort -t, -k2 -n | head -1 | cut -d, -f1 | tr -d '[:space:]')
export CUDA_VISIBLE_DEVICES="${FREE_GPU}"
echo "[tests] using GPU ${CUDA_VISIBLE_DEVICES}"

# Run the solver unit tests (label solver) serially for stable GPU memory.
cd "${BUILD}"
"${CTEST}" -L solver --output-on-failure -j1 2>&1 | tail -60
