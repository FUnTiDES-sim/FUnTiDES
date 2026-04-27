#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status,
# if an unset variable is used, or if a command in a pipeline fails.
set -euo pipefail

# --- Configuration ---
BIN_DIR="./bin/benchmarks"
OUT_DIR="./results"
FORMAT="json"

# Allow passing the number of threads as an argument. Default is 1.
THREADS=${1:-1}

# List of benchmarks to execute
BENCHMARKS=(
    "bench_solver_fe_cartesian_struct_acoustic"
    "bench_solver_fe_cartesian_unstruct_acoustic"
    "bench_solver_fe_cartesian_struct_elastic"
    "bench_solver_fe_cartesian_unstruct_elastic"
)

# --- Execution ---
echo "=========================================================="
echo "Starting the benchmark suite"
echo "Number of Kokkos threads : ${THREADS}"
echo "Output directory        : ${OUT_DIR}/"
echo "=========================================================="

# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

for bench in "${BENCHMARKS[@]}"; do
    echo "[$(date +'%H:%M:%S')] Running: $bench..."

    # Execute and redirect the output to the configured directory
    "$BIN_DIR"/"$bench" \
        --kokkos-threads="$THREADS" \
        --benchmark_format="$FORMAT" \
        > "$OUT_DIR"/"$bench"."$FORMAT"
done

echo "=========================================================="
echo "All benchmarks executed successfully."
echo "=========================================================="
