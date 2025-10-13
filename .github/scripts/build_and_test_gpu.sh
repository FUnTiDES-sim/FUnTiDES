#!/usr/bin/env bash

# Ensure the script fails on errors and undefined variables
#set -euo pipefail

# Expects environment variables already exported:
# SLURM_PARTITION, PYWRAP, PROGRAMMING_MODEL

echo "[build]"
echo "  Partition      = $SLURM_PARTITION"
echo "  Pywrap         = $PYWRAP"
echo "  Model          = $PROGRAMMING_MODEL"

if [[ "$SLURM_PARTITION" == "maple_mig" || "$SLURM_PARTITION" == "maple" ]]; then
  source scripts/environments/env_Maple_GH200.sh
  module list
fi

if [[ "$PYWRAP" == "pywrap-on" ]]; then
  python3 -m venv .venv
  source .venv/bin/activate
  pip install --upgrade pip
  pip install pybind11 pytest numpy
fi

CMAKE_FLAGS=" -DCOMPILE_SEM=ON -DCOMPILE_FD=ON -DENABLE_CUDA=ON"

if [[ "$PYWRAP" == "pywrap-on" ]]; then
  CMAKE_FLAGS+=" -DENABLE_PYWRAP=ON"
else
  CMAKE_FLAGS+=" -DENABLE_PYWRAP=OFF"
fi
[[ "$PROGRAMMING_MODEL" == "vector" ]] && CMAKE_FLAGS+=" -DUSE_VECTOR=ON"
[[ "$PROGRAMMING_MODEL" == "kokkos" ]] && CMAKE_FLAGS+=" -DUSE_KOKKOS=ON"

mkdir -p build
cd build
cmake .. ${CMAKE_FLAGS} -DCMAKE_INSTALL_PREFIX=../install
make -j"$(nproc)"
#ctest --output-on-failure TODO fix on GPU
make install

if [[ "$PYWRAP" == "pywrap-on" ]]; then
  source ../.venv/bin/activate
  export INSTALL_DIR
  INSTALL_DIR=$(realpath ../install)
  export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_DIR/lib
  OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 KOKKOS_NUM_THREADS=2 pytest ../tests/units
fi
