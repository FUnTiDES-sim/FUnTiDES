module load cmake gcc/13.1.1 cuda/12.6.20 

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

export Kokkos_DIR=/shared/data1/Users/j0535952/work2025/ProgrammingModels/kokkos/kokkos/install_maple_rdc/lib64/cmake/Kokkos
export CUDA_ROOT=/hrtc/apps/cuda/12.6.20/aarch64/rocky9/
export CUDA_ARCHITECTURES=90
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1

# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .. ; make
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_OMP=ON  .. ; make

# USE_SEMOPTIM
# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DUSE_SEMOPTIM=ON -DCMAKE_BUILD_TYPE=Release ..; make -j 72 

# USE_SHIVA=ON 
# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DUSE_SHIVA=ON -DCMAKE_BUILD_TYPE=Release ..; make -j 72

# USE_SEMGEOS=ON
# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DUSE_SEMGEOS=ON -DCMAKE_BUILD_TYPE=Release ..;  make -j 72

