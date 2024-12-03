export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
module load cmake/3.26.4 gcc

export DEVICE=p3_Power9CPU
export OMP_PROC_BIND=spread; export OMP_PLACES=threads

# for omp on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make
# for sequential mode on CPU:  cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make
