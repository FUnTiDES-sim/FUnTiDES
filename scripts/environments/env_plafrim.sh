# Usage: FUNTIDES_INSTALL=<FUnTiDES install dir> source env_plafrim.sh

# Clear existing modules to prevent conflicts
module purge

# Load your specific modules
module load tools/git/2.36.0
module load compiler/cuda/12.3
module load language/python/3.12.1
module load compiler/gcc/15.1.0
module load build/cmake/3.27.0
module load mpi/openmpi/5.0.1

# Adding funtides python and lib
if [ -n "${FUNTIDES_INSTALL}" ]; then
  export PYTHONPATH=$FUNTIDES_INSTALL/python${PYTHONPATH:+:$PYTHONPATH}
  export LD_LIBRARY_PATH=$FUNTIDES_INSTALL/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
else
  echo "env_plafrim.sh: FUNTIDES_INSTALL not set, PYTHONPATH and LD_LIBRARY_PATH unchanged" >&2
fi
