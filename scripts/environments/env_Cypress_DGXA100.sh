module av
source "/apps/modules/modulefiles3/supporting_scripts/modver_config.sh"

#Don't use spider cache. This hopefully makes modules v2 and v3 work nicely together
export LMOD_IGNORE_CACHE=1

# activate modules 3
source $SUPPORTSCRIPT_3_DIR/activateMod3
module av
module load cmake gcc/11.4.1 cuda/12.6.20
module list

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
