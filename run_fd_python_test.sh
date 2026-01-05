#!/bin/bash
# Helper script to run FD solver Python tests with proper environment setup

# Set up Python and library paths
export PYTHONPATH=buildPy/python:$PYTHONPATH
export LD_LIBRARY_PATH=buildPy/src/solver/fd:buildPy/src/model/grid:buildPy/src/discretization/fd/kernels:buildPy/src/discretization/fd/stencils:buildPy/src/discretization/fd/abckernels:buildPy/src/utils:$LD_LIBRARY_PATH

# Set up Kokkos/CUDA environment for best performance
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export CUDA_VISIBLE_DEVICES=0
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1

# Run the test
if [ -z "$1" ]; then
    echo "Usage: $0 <test_script>"
    echo "Examples:"
    echo "  $0 tests/benchmarks/python/test_fd_simple.py"
    echo "  $0 tests/benchmarks/python/test_fd_solver_example.py"
    exit 1
fi

python3 "$@"
