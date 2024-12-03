module load cmake rocm/6.2.2

export CC=/opt/rocm/llvm/bin/clang
export CXX=/opt/rocm/llvm/bin/clang++

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
export DEVICE=Mi300aeval_AMDZen4CPU

# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make

