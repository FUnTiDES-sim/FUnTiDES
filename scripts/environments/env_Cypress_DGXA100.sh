module load cmake/3.31.8 gcc/13.2.0 cuda/12.6.20
module list

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
