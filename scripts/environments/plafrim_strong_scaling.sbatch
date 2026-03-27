#!/bin/bash
#SBATCH --job-name=strong-scaling
#SBATCH --output=strong_%j.txt
#SBATCH --error=strong_err_%j.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --nodelist=sirocco[07-13]
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH --exclusive

module purge
module load compiler/cuda/12.3
module load compiler/gcc/12.2.0

export FUNTIDES_PATH=/beegfs/abandet2/projects/FUnTiDES
export FUNTIDES_INSTALL=/beegfs/abandet2/install/funtides
export FUNTIDES_BUILD=$FUNTIDES_PATH/build
export SCRIPTS_FOLDER=/beegfs/abandet2/scripts

source $SCRIPTS_FOLDER/load_cuda_mpi.sh
source $FUNTIDES_PATH/venv/bin/activate

export PYTHONPATH=$FUNTIDES_INSTALL/python:$FUNTIDES_BUILD/external/pykokkos-base/kokkos
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FUNTIDES_INSTALL/lib64
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# Reset CSV
rm -f results_strong.csv

# Strong Scaling Params (Fits in 1 GPU)
GLOBAL_EX=400
EY=300
EZ=300
STEPS=1000

echo "=== STARTING STRONG SCALING ==="

for TASKS in 1 2 4 8; do
    echo "Running with $TASKS GPUs"

    mpirun -np $TASKS \
           --bind-to socket \
           python3 -u $FUNTIDES_PATH/examples/fe/solver_multigpu_multinode.py \
           --ex $GLOBAL_EX --ey $EY --ez $EZ --steps $STEPS \
           --benchmark --type-scaling strong

done
