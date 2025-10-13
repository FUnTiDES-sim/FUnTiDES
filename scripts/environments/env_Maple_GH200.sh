source "/apps/modules/modulefiles3/supporting_scripts/modver_config.sh"

#Don't use spider cache. This hopefully makes modules v2 and v3 work nicely together
export LMOD_IGNORE_CACHE=1

# activate moddules 3
source $SUPPORTSCRIPT_3_DIR/activateMod3
module av

module load cmake gcc/13.1.1 cuda/12.6.20
module list

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

export Kokkos_DIR=/shared/data1/Users/j0535952/work2025/ProgrammingModels/kokkos/kokkos/install_maple_rdc/lib64/cmake/Kokkos
export CUDA_ROOT=/hrtc/apps/cuda/12.6.20/aarch64/rocky9/
export CUDA_ARCHITECTURES=90
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1

# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .. ; make