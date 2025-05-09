module load cmake gcc/13.1.1 cuda/12.6.20 

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

export Kokkos_DIR=/shared/data1/Users/j0535952/work2025/ProgrammingModels/kokkos/kokkos/install_maple_rdc/lib64/cmake/Kokkos
export CUDA_ROOT=/hrtc/apps/cuda/12.6.20/aarch64/rocky9/

# by default: USE_SEMOPTIM is ON
#
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .. ; make
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_OMP=ON  .. ; make
# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make 


# set -DUSE_SEMCLASSIC=ON -DUSE_SEMOPTIM=OFF
# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DUSE_SEMCLASSIC=ON -DUSE_SEMOPTIM=OFF .. ; make 

# set -DUSE_SHIVA=ON -DUSE_SEMOPTIM=OFF
# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DUSE_SHIVA=ON -DUSE_SEMOPTIM=OFF .. ; make 


