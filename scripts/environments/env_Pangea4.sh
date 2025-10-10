#!/bin/bash

# use lastest GNU compiler
module purge
module load PrgEnv-gnu python/3.11.6 cmake/3.27.2 craype-x86-milan
module unload cray-libsci cray-mpich
module list