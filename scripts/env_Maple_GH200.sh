module load cmake gcc/13.1.1 cuda/12.6.20 

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

export Kokkos_DIR=/shared/data1/Users/j0535952/work2025/ProgrammingModels/kokkos/kokkos/install_maple_rdc/lib64/cmake/Kokkos
export CUDA_ROOT=/hrtc/apps/cuda/12.6.20/aarch64/rocky9/
export CUDA_ARCHITECTURES=90
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1

# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .. ; make
# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -Duse OMP=ON  .. ; make

# for kokkos on GPU: cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -Duse KOKKOS=ON -DENABLE_CUDA=ON -DCMAKE_BUILD_TYPE=Release ..; make -j 72 

# for use SEMCLASSIC
# ./bin/semproxy  -physics 0  -method 0 -order 2

# for use SEMOPTIM
# ./bin/semproxy  -physics 0  -method 1 -order 2

# for use SEMGEOS
# ./bin/semproxy  -physics 0  -method 2 -order 2

# for use SHIVA
# ./bin/semproxy  -physics 0  -method 3 -order 2

