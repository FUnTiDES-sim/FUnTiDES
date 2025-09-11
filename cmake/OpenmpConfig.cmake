#-------------------------------------------------------------------
# OpenMP configuration
#-------------------------------------------------------------------

# the solver is using OpenMP+Vector on CPU
if(NOT USE_OMP OR USE_KOKKOS)
  message("-- USE_OMP flag is set to be false or Kokkos activated. OMP is not enabled")
elseif()
  message(STATUS "BUILDING SOLVER including OpenMP+Vector on CPU")
  find_package(OpenMP REQUIRED)
  set(ENABLE_CUDA OFF CACHE STRING "" FORCE)
  set(USE_VECTOR ON CACHE BOOL "" FORCE)
  message("-- USE_OMP flag is set: USE_VECTOR= ${USE_VECTOR} ")
endif()