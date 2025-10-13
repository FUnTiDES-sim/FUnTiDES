module av
module load cmake gcc/11.4.1 cuda/12.6.20
module list

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
