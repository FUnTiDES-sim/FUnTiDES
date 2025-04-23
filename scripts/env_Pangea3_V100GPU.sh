export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
module load cmake/3.26.4 gcc cuda

export SEM_TPL_ROOT_DIR=/appli_RD/JIEMENG/SEMproxy/proxy_20240430/tpl4ProxyApp/installTPL_P3_GPU/
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include
export CUDA_ROOT=/data_local/sw/cuda/11.0.3-rhel8/

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

# don't forget to enable "-DPower9_pangea3=ON"

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DPower9_pangea3=ON .. ; make 
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON -DPower9_pangea3=ON .. ; make 
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON -DPower9_pangea3=ON .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DPower9_pangea3=ON .. ; make

