module load cmake gcc/13.1.1 cuda/12.6.20 

export OMP_PROC_BIND=spread; export OMP_PLACES=threads

export CUDA_ROOT=/hrtc/apps/cuda/12.6.20/aarch64/rocky9/
export CUDA_ARCHITECTURES=90
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1

# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .. ; make

# for kokkos on GPU:
# Optional MPS setup to avoid conflicts from other users
# export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-$USER
# export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-mps-log-$USER
# mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
# nvidia-cuda-mps-control -d
# Then run cmake command with target architecture and build type, e.g.:
# cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON -DCMAKE_BUILD_TYPE=Release ..; make -j 72 
