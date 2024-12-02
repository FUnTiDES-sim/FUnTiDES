cmake_minimum_required( VERSION 3.22.1 )

project( proxyapp )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------

set (CMAKE_VERBOSE_MAKEFILE ON )

#-----------------------------------------------------
# Setup the computing resources for running proxyapp
#-----------------------------------------------------
option (X86_laptop "X86 architecture  default" ON)
option (X86_A100 "Compilation on X86CPU+A100 node" OFF)
option (Power9_pangea3 "Power9+V100 node" OFF)
option (ARM "ARM aarch64 architecture: e.g. grace CPU" OFF)

#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------

# use OMP to parallelise loops
option (USE_OMP "Use OpenMP to parallelise loops" OFF)

# use RAJA to expose parallelism loop
option (USE_RAJA "Use RAJA to parallelise loops" OFF)

# use KOKOS to expose parallelism loop
option (USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)

option (USE_SEMCLASSIC "Use CLASSIC " OFF)
option (USE_SEMOPTIM "Use OPTIM " OFF)
option (USE_SHIVA "Use SHIVA " OFF)

#---------------------------------------------------
# Tests directory to perform differen tests before
# integration into proxy apps
#---------------------------------------------------
option (addTests "add Tests directory" OFF)

#--------------------------------------------------------
# Option to save snapshots of wavefields (for debugging)
#--------------------------------------------------------
option (SEM_SAVE_SNAPSHOTS "save snapshots for SEM-proxy" OFF)
option (FD_SAVE_SNAPSHOTS "save snapshots for FD-proxy" OFF)

#------------------------------------------------------
# Option to printout memory allocation info of arrays
#------------------------------------------------------
option (PRINT_ALLOC_INFO "printout memory allocation info" OFF)

#---------------------------------------------------------
# Option to use hierarchical parallelism in Kokkos
#---------------------------------------------------------
option (USE_KOKKOS_TEAMS "use hierarchical parallelism in Kokkos" OFF)

add_subdirectory(src)


