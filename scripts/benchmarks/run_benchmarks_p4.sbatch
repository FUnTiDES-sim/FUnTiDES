#!/bin/sh
  
#SBATCH --job-name="benchmark_funtides"
#SBATCH --comment=MAKUTU,SSC,benchmarks
#SBATCH --time=2:00:00
#SBATCH --partition=p4_dev
#SBATCH --output=benchmark_%j.out
#SBATCH --error=benchmark_%j.err
##SBATCH --exclusive
#SBATCH --mem=734GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1 # avoid hyperthreading
#SBATCH --sockets-per-node=1 # allocate one socket
#SBATCH --cpus-per-task=96   # allocate the full socket

ulimit -s unlimited
ulimit -c unlimited

# SETUP
# assumes
#  - script is in scripts/benchmarks and submitted from repo root
#  - build & install dirs are in repo root
#  - environment script is in scripts/environments
#  - virtual env is in repo root .venv
ROOT_DIR=$SLURM_SUBMIT_DIR
SCRIPT_DIR=$ROOT_DIR/scripts/benchmarks
ENV_DIR=$ROOT_DIR/scripts/environments
BUILD_DIR="$ROOT_DIR/build"
INSTALL_DIR="$ROOT_DIR/install"

echo "=== SETUP SCRIPT ==="
echo "Script dir: $SCRIPT_DIR"
echo "Root dir: $ROOT_DIR"
echo "Build dir: $BUILD_DIR"
echo "Install dir: $INSTALL_DIR"

source $ENV_DIR/env_Pangea4.sh

source $ROOT_DIR/.venv/bin/activate
export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_DIR/lib64

# Disclaimer -> OpenMP affinity setup is probably not optimal yet, just an attempt
export OMP_PLACES=cores
export OMP_PROC_BIND=close # close (cache locality) | spread (bandwidth / NUMA balance)
export OMP_DYNAMIC=FALSE

echo "=== SETUP OMP ==="
echo "OMP_PLACES=$OMP_PLACES"
echo "OMP_PROC_BIND=$OMP_PROC_BIND"
echo "OMP_DYNAMIC=$OMP_DYNAMIC"

# CHECK
echo "=== TOPOLOGY ==="
lscpu | egrep 'Socket|Core|CPU\(s\)|NUMA'
echo
echo "numactl -H:"
numactl -H
echo "=== BINDING ==="
echo
echo "SLURM env:"
echo "SLURM_JOB_CPUS_PER_NODE=$SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE"
echo
srun --cpu-bind=socket,verbose --cpus-per-task=96 bash -c '
echo "Allocated mask:"; grep Cpus_allowed_list /proc/self/status
echo "Taskset:"; taskset -cp $$
echo "numactl --show:"; numactl --show
'

# RUN
echo "=== PYTHON BENCHMARKS ==="
pushd $ROOT_DIR
srun --cpu-bind=socket,verbose --cpus-per-task=96 python scripts/run_pywrap_benchmarks.py --threads 2,4,8,16,32,64
popd
echo "=== CPP BENCHMARKS ==="
pushd $BUILD_DIR
srun --cpu-bind=socket,verbose --cpus-per-task=96 ctest -L benchmark
popd