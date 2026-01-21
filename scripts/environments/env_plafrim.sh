# Clear existing modules to prevent conflicts
module purge

# Load your specific modules
module load tools/git/2.36.0
module load compiler/cuda/12.3
module load language/python/3.12.1
module load compiler/gcc/15.1.0
module load build/cmake/3.27.0
module load mpi/openmpi/5.0.1
