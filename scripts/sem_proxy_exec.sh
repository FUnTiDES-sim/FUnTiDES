#!/bin/bash
# Run wrapper for the SEM proxy harnesses, placing the benchmark on the
# compute node.
#
# MIMIR runs on the Oleron login node (glibc 2.28), which cannot load the
# FUnTiDES stack (GLIBC 2.34 required). elba-b02 (partition COMPUTE_GPU_B300,
# RHEL 9, glibc 2.34, 8x B300) can. This wrapper is registered as the proxy
# executable_path, so every invocation the server makes locally --
# proxy_exec(op='reference') and proxy runs -- is executed on elba-b02
# through srun.
#
# Location-aware: inside a Slurm job that already placed us on the node (the
# proxy_slurm route), exec the binary directly instead of nesting a second
# srun inside the allocation.
#
# Usage: sem_proxy_exec.sh <acoustic|elastic> <output-file> [extra args...]
#
# MEASUREMENT HYGIENE -- why this wrapper waits. elba-b02 is a shared node:
# neighbours without a gres allocation have been observed holding GPUs. A
# contended sample inflates per-step time by 3x-16x (measured 2026-09), which
# would make a refactor's no-regression gate accept noise. So the run waits
# until the selected device is quiet. This changes WHEN the benchmark is
# measured, never WHAT it measures: the binary, its flags and its timer are
# untouched.
#
# GPU selection: prefer the device Slurm allocated (CUDA_VISIBLE_DEVICES set
# by the job); otherwise pick the device with the least memory in use. The
# harness initialises Kokkos without a device argument, so the first visible
# device is the one used.
#
# Overridable without editing this file:
#   SEMPROXY_NODE          node to run on              (default elba-b02)
#   SEMPROXY_PARTITION     partition that selects it   (default COMPUTE_GPU_B300)
#   SEMPROXY_QUIET_UTIL    utilisation ceiling, percent (default 5)
#   SEMPROXY_QUIET_WAIT    maximum seconds to wait      (default 45)
#   SEMPROXY_QUIET_ON       set 0 to disable the wait   (default 1)
#
# The wait default is 45s, not 900s: the proxy server's blocking reference
# seal must complete within its 120s tool budget (configure ~25s + case ~6s),
# and a longer wait killed the first seal attempt mid-flight (the case ran
# but the field was never sealed). If the GPU is still busy after the wait,
# the run proceeds with a warning: l2_rel (the hard gate) is deterministic
# and contention-immune, time_s is interpreted with the machine noise floor.
set -uo pipefail

PHYSICS=${1:?usage: sem_proxy_exec.sh <acoustic|elastic> <output-file> [extra args]}
shift

NODE=${SEMPROXY_NODE:-elba-b02}
PARTITION=${SEMPROXY_PARTITION:-COMPUTE_GPU_B300}
QUIET_UTIL=${SEMPROXY_QUIET_UTIL:-5}
QUIET_WAIT=${SEMPROXY_QUIET_WAIT:-45}
QUIET_ON=${SEMPROXY_QUIET_ON:-1}

case "$PHYSICS" in
  acoustic)     BINARY=proxy_sem_acoustic ;;
  elastic)      BINARY=proxy_sem_elastic ;;
  elastic_tti)  BINARY=proxy_sem_elastic_tti ;;
  *) echo "[exec] unknown physics '$PHYSICS' (expected acoustic|elastic|elastic_tti)" >&2; exit 2 ;;
esac
REAL=/appli_RD/j1039830/FUnTiDES/build/bin/benchmarks/${BINARY}

# Force the Kokkos this binary was built against (sm_103) onto the loader
# path. Both Kokkos installs export the same SONAME (libkokkoscore.so.5.1) and
# the loader consults LD_LIBRARY_PATH before DT_RUNPATH, so pinning makes the
# outcome deterministic. Pairs with SEMPROXY_KOKKOS_DIR in build_sem_proxy.sh.
KOKKOS_PREFIX=${SEMPROXY_KOKKOS_PREFIX:-/appli_RD/j1039830/local}
export LD_LIBRARY_PATH="${KOKKOS_PREFIX}/lib64:${LD_LIBRARY_PATH:-}"

gpu_util() {
  nvidia-smi -i "$1" --query-gpu=utilization.gpu --format=csv,noheader,nounits \
    2>/dev/null | tr -d '[:space:]'
}

pick_gpu() {
  # Prefer the device Slurm allocated; fall back to the least-used one.
  local slum_gpu=${CUDA_VISIBLE_DEVICES:-}
  if [ -n "$slum_gpu" ] && [ "${slum_gpu#*,}" = "$slum_gpu" ]; then
    echo "$slum_gpu"
    return 0
  fi
  nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
    | sort -t, -k2 -n | head -1 | cut -d, -f1 | tr -d '[:space:]'
}

wait_for_quiet_gpu() {
  [ "$QUIET_ON" = "1" ] || return 0
  local gpu="$1"
  local waited=0 util
  while [ "$waited" -lt "$QUIET_WAIT" ]; do
    util=$(gpu_util "$gpu")
    if [ -z "$util" ]; then
      echo "[exec] WARNING: cannot read GPU $gpu utilisation; running anyway" >&2
      return 0
    fi
    if [ "$util" -lt "$QUIET_UTIL" ]; then
      # One confirming sample: the neighbour may be between kernels.
      sleep 2
      util=$(gpu_util "$gpu")
      if [ -n "$util" ] && [ "$util" -lt "$QUIET_UTIL" ]; then
        echo "[exec] GPU $gpu quiet at ${util}% after ${waited}s wait" >&2
        return 0
      fi
    fi
    sleep 10
    waited=$((waited + 10))
  done
  echo "[exec] WARNING: GPU $gpu still at ${util}% after ${QUIET_WAIT}s;" \
       "running anyway -- this sample may be contended" >&2
  return 0
}

if [ "${SEMPROXY_ON_NODE:-0}" = "1" ] || [ "$(hostname -s)" = "$NODE" ]; then
  GPU=$(pick_gpu)
  export CUDA_VISIBLE_DEVICES="$GPU"
  wait_for_quiet_gpu "$GPU"
  exec "$REAL" "$@"
fi

exec /usr/bin/srun --partition="$PARTITION" --nodelist="$NODE" --gres=gpu:1 \
  --cpus-per-task=8 --mem=32G --time=01:00:00 --quiet \
  --comment="FUnTiDES,LLM4CSE,OT/R&D/SUB/SC" \
  env SEMPROXY_ON_NODE=1 \
  /appli_RD/j1039830/FUnTiDES/scripts/sem_proxy_exec.sh "$PHYSICS" "$@"
