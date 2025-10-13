module av
export MODULEPATH=/apps/modules/modulefiles3/Linux:/apps/modules/modulefiles3/x86_64/DGX/Common:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/Tools:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/SciLib:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/Debugger:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/MPI:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/Compilers:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/Compilers_spack:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/Sched:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/PrgEnv:/apps/modules/modulefiles3/x86_64/DGX/Core/Production/Common:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/Tools:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/HPC/MPI:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/HPC/Tools:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/HPC/Compiler:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/SciLib:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/Debugger:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/MPI:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/Compilers:/apps/modules/modulefiles3/x86_64/DGX/Core/Experimental/Sched:/apps/modules/modulefiles3/Common/proxy:/home/svcmakutu/user_env_files
#Don't use spider cache. This hopefully makes modules v2 and v3 work nicely together
export LMOD_IGNORE_CACHE=1

module av
module load cmake gcc/11.4.1 cuda/12.6.20
module list

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
