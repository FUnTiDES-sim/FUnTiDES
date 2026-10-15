#!/bin/bash
# Build wrapper for the SEM proxy harnesses (proxy_sem_acoustic, proxy_sem_elastic).
#
# The proxy server runs build_cmd as argv (not a shell line), so the PATH setup
# a shell would normally provide lives here.
#
# BUILD WITHOUT -O3: the SEM solver has a latent use-after-free that -O3
# (CMAKE_BUILD_TYPE=Release) turns into a crash in
# SEMsolver::updateFieldsForward (Kokkos ReferenceCountedDataHandle copy).
# The reference state is the non-O3 build. Do not add -O3.
set -euo pipefail

# MIMIR runs on the Oleron login node (glibc 2.28), which cannot run the
# project's Python 3.11 venv -- `cmake` configure therefore dies in
# find_package(Python REQUIRED) there. Build on the compute node instead: if
# we are not on the node, re-invoke this same script on it through srun.
# Inside a Slurm job that already placed us on the node (the proxy_slurm
# route) this is a no-op and the build proceeds directly.
#
# The build MUST run on the node the binary will execute on: CMAKE_CUDA_ARCHITECTURES
# is "native" in the project, so the CUDA target is detected from the GPU
# present at configure time. elba-b02 is the B300 node (sm_103).
#
# Overridable without editing this file:
#   SEMPROXY_BUILD_NODE      node to build on            (default elba-b02)
#   SEMPROXY_BUILD_PARTITION partition that selects it    (default COMPUTE_GPU_B300)
#   SEMPROXY_KOKKOS_DIR      Kokkos CMake dir            (default local, sm_103)
NODE=${SEMPROXY_BUILD_NODE:-elba-b02}
PARTITION=${SEMPROXY_BUILD_PARTITION:-COMPUTE_GPU_B300}
if [ "$(hostname -s)" != "$NODE" ]; then
  exec /usr/bin/srun --partition="$PARTITION" --nodelist="$NODE" \
    --cpus-per-task=32 --mem=64G --time=02:00:00 --quiet \
    --comment="mimir proxy build SEM harnesses (FUnTiDES solver refactor)" \
    /appli_RD/j1039830/FUnTiDES/scripts/build_sem_proxy.sh
fi

export PATH=/appli_RD/j1039830/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/appli_RD/j1039830/local/cuda/lib64:${LD_LIBRARY_PATH:-}

# Configure on the node, where FindPython can execute the venv interpreter.
# Idempotent: a few seconds when the cache and the CMake inputs are unchanged.
#
# Kokkos_DIR is pinned explicitly. The installed Kokkos is an arch LOCK, not a
# default: its KokkosConfigCommon.cmake sets Kokkos_CUDA_ARCHITECTURES and
# propagates it to every consumer, so this is what selects the target GPU:
#   sm_103 -> /appli_RD/j1039830/local/lib64/cmake/Kokkos          (B300, elba-b02)
#   sm_120 -> /appli_RD/j1039830/local/tpl-sm120/lib64/cmake/Kokkos (RTX 6000, elba-x02)
#
# CMAKE_BUILD_TYPE is forced EMPTY: a cached Release from an earlier session
# would silently add -O3 and the build would produce a crashing binary.
CMAKE=/appli_RD/j1039830/local/cmake/bin/cmake
SRC=/appli_RD/j1039830/FUnTiDES
BUILD=/appli_RD/j1039830/FUnTiDES/build
KOKKOS_DIR=${SEMPROXY_KOKKOS_DIR:-/appli_RD/j1039830/local/lib64/cmake/Kokkos}
"${CMAKE}" -S "${SRC}" -B "${BUILD}" \
  -DCMAKE_BUILD_TYPE= \
  -DPython_EXECUTABLE=/appli_RD/j1039830/local/tpl-venv/bin/python \
  -DKokkos_DIR="${KOKKOS_DIR}"

# Parallel build: fun_solver is a dozen heavy template TUs that rebuild on
# every header touch; serial make costs ~20 min. 192 cores are available; cap
# at 32 to keep memory per compiler sane.
make -j32 -C "${BUILD}" proxy_sem_acoustic proxy_sem_elastic proxy_sem_elastic_tti
