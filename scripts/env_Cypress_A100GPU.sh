module load gcc/12.2.0 cuda/12.0.76  cmake

export SEM_TPL_ROOT_DIR=/shared/data1/Users/j0535952/work2024/proxys_2024/tpl4ProxyApp/installTPL_Cypress_A100GPU

export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include
export CUDA_ROOT=/hrtc/apps/cuda/12.0.76/x86_64/centos8

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make install
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make install
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make

