module load cmake gcc/13.1.1 cuda/12.6.20
module list

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

export Kokkos_DIR=/shared/data1/Users/j0535952/work2025/ProgrammingModels/kokkos/kokkos/install_maple_rdc/lib64/cmake/Kokkos
export CUDA_ROOT=/hrtc/apps/cuda/12.6.20/aarch64/rocky9/
export CUDA_ARCHITECTURES=90
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1