## Welcome to the  proxyApp project!

The proxyApp project  collects a suite of simple codes representing real applications.
It is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.


## Actual applications 

Current implementation of the proxyApp includes SEM ( Spectral finite Element Methods) and FD ( Finite Differences methods) to solve 2nd order acoustic wave equation in 2D and 3D spaces:
* The SEM proxy applicaton is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations.
* The FD proxy applicaton is a benchmark designed to simulate wave propagation using finite differences stecnils operators for solving partial differential equations.

One of the key features of the SEM and FD proxy benchmarks are their adaptability to different programming models and HPC architectures. This makes them a useful proxy applications for advancing the state of the art in high-performance computing. In addition to their technical capabilities, they are designed to be easy to build and use. This makes them accessible to a wide range of users, from researchers to developers.

## What Programming Models users can select to run  proxy applications?

The programming models included in the current sem proxy implementations include:
* OpenMP [https://www.openmp.org/] to parallelizing for loops
* RAJA [https://raja.readthedocs.io/en/develop/]
* KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]

RAJA and Kokkos must first be compiled and installed  in TPL4ProxyApp repository.

## What data containers users can select to run SEM proxy?

The data containers included in the current sem proxy implementations include:
* LvArray [https://lvarray.readthedocs.io/en/latest/]
* C++ std::vector

## Quick Start to compile and install:

### Step 1: compile and install proxyApp

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install  ..  
   make; make install
```

The default compilation is sequential mode with std::vector implementation. 
So you will get an executable named "sem_SEQUENTIAL.exe" and FDTDSEQUENTIAL.exe in your installation directory.

### Step 2: run the executable, for example:

```
   install/bin/proxyName_SEQUENTIAL.exe ( with proxyName: sem or fd)
```

## Available configuration:

### OPEN_MP
in the case of OPENMP std::vector container is used.
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON  ..  
   make; make install
```
### RAJA + OPEN_MP + CUDA (on Nvidia GPUs)
in the case of RAJA Lvarray container is used.
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON  ..  
   make; make install
 
```

### KOKKOS + OPEN_MP + CUDA (on Nvidia GPUs)
KOKKOS provides its own data container.
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON  ..  
   make install
```

### RAJA + OPEN_MP + HIP (on AMD GPUs)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON  ..  
   make; make install
```

### KOKKOS + OPEN_MP + HIP (on AMD GPUs)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON  ..  
   make install
 
```

### RAJA + OPEN_MP + CUDA + ARM (on Nvidia Grace-Hopper)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON ..  
   make; make install
```

### KOKKOS + OPEN_MP + CUDA + ARM (on Nvidia Grace-Hopper)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON  -DARM=ON ..  
   make install
```
