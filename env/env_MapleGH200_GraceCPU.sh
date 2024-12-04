module load cmake gcc/13.1.1 cuda

export CC=/opt/rh/gcc-toolset-13/root/bin/gcc
export CXX=/opt/rh/gcc-toolset-13/root/bin/g++

export DEVICE=maple_GraceCPU

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make

