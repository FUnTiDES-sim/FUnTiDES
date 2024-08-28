#!/bin/sh

# to compile 

source env/env_{system}.sh

mkdir build; cd build

# sem
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DCMAKE_BUILD_TYPE=RelWithDebInfo .. ; make 

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make sem_Raja.exe
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make sem_Kokkos.exe

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON .. ; make sem_Raja.exe
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON .. ; make sem_Kokkos.exe 

# fd

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make fd_Raja.exe
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make fd_Kokkos.exe

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON .. ; make fd_Raja.exe
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON .. ; make fd_Kokkos.exe

# on grace-hopper

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DARM=ON .. ; make 

# on pangea3

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON -DPower9_pangea3=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DPower9_pangea3=ON .. ; make 
