module load cmake rocm/6.2.2

export CC=/opt/rocm/llvm/bin/clang
export CXX=/opt/rocm/llvm/bin/clang++

export KOKKOS_DEV=/shared/data1/Users/j0535952/work2024/ProgrammingModels/kokkos_Dec2024/kokkos/
export KOKKOS_DIR=${KOKKOS_DEV}/build_MI300A_clang++/install/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${KOKKOS_DEV}/include

export DEVICE=Mi300aeval_Mi300GPU
export OMP_PROC_BIND=spread; export OMP_PLACES=threads
export HSA_XNACK=1;

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON .. ; make 
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON .. ; make 
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make

