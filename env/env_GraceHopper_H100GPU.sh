module load gcc/12.2.0 cmake cuda/12.4.131

export SEM_TPL_ROOT_DIR=/shared/data1/Users/j0535952/work2024/proxys_2024/tpl4ProxyApp_essential/tpl4ProxyApp/installTPL_Maple_GPU
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/chai
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib64/cmake/umpire/
export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include
export CUDA_ROOT=/hrtc/apps/cuda/12.4.131/aarch64/rocky9/

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DARM=ON -DX86_laptop=OFF .. ; make 
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON -DX86_laptop=OFF .. ; make 
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON -DARM=ON -DX86_laptop=OFF .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DARM=ON -DX86_laptop=OFF .. ; make
