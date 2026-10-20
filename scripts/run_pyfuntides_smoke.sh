#!/bin/bash
# Smoke test: pyfuntides import + one GPU solver run through the Python API.
# Validates that the fe/api boundary and the pywrap bindings survive the
# src/solver refactor. Runs on elba-b02 (glibc 2.34 + GPU required).
#
# Known pre-existing issues deliberately OUT of scope (see workspace memory):
# - n_rhs=1 crashes (RHSElement[1] OOB) -> we use n_rhs=2
# - tests/validation conftest --debug conflicts with pytest 9.1.1
# - tests/integration hardcode HostSpace/LayoutRight (CUDA expects UVM/Left)
set -euo pipefail

NODE=${SEMPROXY_BUILD_NODE:-elba-b02}
PARTITION=${SEMPROXY_BUILD_PARTITION:-COMPUTE_GPU_B300}
if [ "$(hostname -s)" != "$NODE" ]; then
  exec /usr/bin/srun --partition="$PARTITION" --nodelist="$NODE" \
    --gres=gpu:1 --cpus-per-task=16 --mem=48G --time=02:00:00 --quiet \
    --comment="FUnTiDES,LLM4CSE,OT/R&D/SUB/SC" \
    /appli_RD/j1039830/FUnTiDES/scripts/run_pyfuntides_smoke.sh
fi

source /appli_RD/j1039830/local/tpl-venv/bin/activate
export PATH=/appli_RD/j1039830/local/cmake/bin:/appli_RD/j1039830/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/appli_RD/j1039830/local/cuda/lib64:/appli_RD/j1039830/local/lib64:${LD_LIBRARY_PATH:-}

# Rebuild the pywrap module (links the refactored libfun_solver) and copy the
# freshly built module into the installed package so the smoke exercises the
# CURRENT tree. A full `cmake --install` is not used: this build tree has never
# built every installable target (e.g. libfun_io.so), which a full install
# requires; the smoke only needs the solver module refreshed.
cd /appli_RD/j1039830/FUnTiDES
cmake --build build --target solver -j16
cp build/src/solver/sem/pywrap/solver.cpython-311-x86_64-linux-gnu.so \
   /appli_RD/j1039830/local/funtides-install/python/pyfuntides/solver.cpython-311-x86_64-linux-gnu.so

export PYTHONPATH="/appli_RD/j1039830/local/funtides-install/python:${PYTHONPATH:-}"

cd /appli_RD/j1039830/FUnTiDES

echo "[smoke] 1/2 import pyfuntides"
python -c "import pyfuntides; import pyfuntides.solver as S; print('import OK:', S.__name__)"

echo "[smoke] 2/2 GPU solver run (solver_cartesian, MAKUTU, n_rhs=2)"
SMOKE_LOG=/appli_RD/j1039830/FUnTiDES/proxy_bench/smoke/solver_cartesian_run.log
python -u examples/fe/solver_cartesian.py --mem GPU --impl MAKUTU --order 1 \
  --ex 8 --ey 8 --ez 8 --n_time_steps 20 --n_rhs 2 > "${SMOKE_LOG}" 2>&1
RUN_RC=$?
echo "[smoke] run exit code: ${RUN_RC} (full log: ${SMOKE_LOG})"
tail -25 "${SMOKE_LOG}"
