module load cmake gcc/13.1.1 cuda

export CC=/opt/rh/gcc-toolset-13/root/bin/gcc
export CXX=/opt/rh/gcc-toolset-13/root/bin/g++

export KOKKOS_DEV=/shared/data1/Users/j0535952/work2024/ProgrammingModels/kokkos_Dec2024/kokkos/
export KOKKOS_DIR=${KOKKOS_DEV}/build_H100_g++_rdc/install/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${KOKKOS_DEV}/include

export CUDA_ROOT=/hrtc/apps/cuda/12.4.131/aarch64/rocky9
export CUDA_ARCHITECTURES=90
export DEVICE=maple_H100GPU

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make 
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make 

