#-------------------------------------------------------------------
# The Solver Implementation is with Programming Models:
#
#  Users can choose from "OMP", "RAJA", and "KOKKOS" as following:
#
#  (0) The default is using sequential mode
#  (1) If USE_OMP is ON, it is using OpenMP+Vector on CPU.
#  (2) If USE_RAJA is ON, it is using RAJA+LVARRAY+CUDA.
#  (3) If USE_KOKKOS is ON, it is using KOKKOS+CUDA.
#
#  And the users can chose data containers as following:
#
#  env variables to be defined before build
#  the variables above are defined for Pangea3
#  _TPL must be first compiled and installed
#
#  example:
#  export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
#  module load cmake/3.21.4 gcc/8.4.1 cuda/11.0.3 ompi/4.1.2 openblas/0.3.18 python4geosx/3.8.5-gcc-8.4.1-ompi-4.1.2
#  export _TPL_ROOT_DIR=/users/j0024947/appli_src/geosx/codes/proxyApp/tpl4ProxyApp/installTPL
#  export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
#  export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
#  export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
#  export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
#
#-------------------------------------------------------------------

# the solver is using sequential mode by default
set (SOLVER "SEQUENTIAL" CACHE PATH "" FORCE)
set (USE_VECTOR ON CACHE BOOL "" FORCE)

if(ENABLE_HIP)
  set(COMP_VER "6.0.0")
  set(COMP_ARCH "gfx908")
  #set(CMAKE_CXX_COMPILER "/opt/rocm-6.0.0/bin/amdclang++")
  set(CMAKE_CXX_COMPILER "/opt/rocm-6.0.0/bin/hipcc")
  set(ROCM_ROOT_DIR "/opt/rocm-${COMP_VER}")
  set(ROCM_PATH "/opt/rocm-${COMP_VER}")
  set(HIP_ROOT_DIR "/opt/rocm-${COMP_VER}/hip")
  set(HIP_PATH "/opt/rocm-${COMP_VER}/bin")
  set(GPU_TARGETS "${COMP_ARCH}")
  message(STATUS "======== ROCM_ROOT_DIR=${ROCM_ROOT_DIR}")
  set(USE_SEM_INLINE ON CACHE BOOL "" FORCE) #clang doesn't support -fgpu-rdc yet
else()
  # CC options
  set(CMAKE_C_COMPILER  /usr/bin/gcc-12  CACHE PATH "")
  set(CMAKE_C_FLAGS_RELEASE "-g -O3 -DNDEBUG" CACHE STRING "")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

  # C++ options
  set(CMAKE_CXX_COMPILER /usr/bin/g++-12  CACHE PATH "")
  message( STATUS, "-----------------------CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
  set (CMAKE_CXX_FLAGS "-g -O3 -std=c++20 -fconcepts " CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -DNDEBUG " CACHE STRING "")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
  set(CMAKE_CXX_STANDARD 20 CACHE STRING "")
endif()

# allow to install only one target
set (CMAKE_SKIP_INSTALL_ALL_DEPENDENCY "True" CACHE STRING "" FORCE)

################################
# BLT
################################

# Default to C++20 if not set so GTest/GMock can build
if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++20" CACHE STRING "")
endif()

# the following options are setup for BLT
set( ENABLE_GTEST OFF CACHE STRING "" FORCE)
set( ENABLE_ALL_WARNINGS OFF CACHE STRING "" FORCE)
set( ENABLE_TESTS OFF CACHE STRING "" FORCE)
set( ENABLE_EXAMPLES OFF )
set( ENABLE_DOCS OFF )
set( ENABLE_BENCHMARKS OFF )
set( DISABLE_UNIT_TESTS ON)
# force ENABLE_OPENMP to ON
set( ENABLE_OPENMP ON CACHE STRING "" FORCE)

# choose platform for Cuda options
if(X86_laptop)
  set(X86_cypress OFF CACHE STRING "" FORCE)
  if (ENABLE_CUDA)
    set(CUDA_ROOT $ENV{CUDA_ROOT} CACHE PATH "" FORCE)
    set(CMAKE_CUDA_ARCHITECTURES 89 CACHE STRING "" FORCE)
    set(CUDA_ARCH sm_89 CACHE STRING "" FORCE)
    #set(CMAKE_CUDA_FLAGS "-lineinfo -g -pg -restrict --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -std=c++20 -fconcepts" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS "-O3 -restrict --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -std=c++20 -fconcepts" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELEASE "-g -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")
  endif()

elseif(Power9_pangea3)
  # Cuda options for pangea3
  if (ENABLE_CUDA)
     set(CUDA_ROOT "/data_local/sw/cuda/11.0.3-rhel8" CACHE PATH "")
     set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "" FORCE)
     set(CUDA_ARCH sm_70 CACHE STRING "" FORCE)
     set(CMAKE_CUDA_FLAGS "-O3 -restrict -arch ${CUDA_ARCH}  --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler -std=c++20 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le  " CACHE STRING "")
     set(DEVICE "V100GPU" CACHE PATH "" FORCE)
  else()
    set(DEVICE "Power9" CACHE PATH "" FORCE)
  endif()

elseif(ARM)
  if (ENABLE_CUDA)
    message(STATUS $ENV{CUDA_ROOT})
    set(CUDA_ROOT $ENV{CUDA_ROOT} CACHE PATH "" FORCE)
    set(CMAKE_CUDA_ARCHITECTURES 90 CACHE STRING "")
    set(CUDA_ARCH sm_90 CACHE STRING "")
    set(CMAKE_CUDA_FLAGS "-restrict --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -std=c++20" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELEASE "-Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")
    set(DEVICE "H100GPU" CACHE PATH "" FORCE)
  else()
    set(DEVICE "GRACE" CACHE PATH "" FORCE)
  endif()

elseif(X86_cypress)
  if (ENABLE_CUDA)
    set(CUDA_ROOT $ENV{CUDA_ROOT} CACHE PATH "" FORCE)
    set(CMAKE_CUDA_ARCHITECTURES 80 CACHE STRING "" FORCE)
    set(CUDA_ARCH sm_80 CACHE STRING "" FORCE)
    #set(CMAKE_CUDA_FLAGS "-O3 -restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler -std=c++20" CACHE STRING "" FORCE)
    set(CMAKE_CUDA_FLAGS "-restrict --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -std=c++20" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELEASE "-Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")
    set(DEVICE "A100GPU" CACHE PATH "" FORCE)
  elseif(ENABLE_HIP)
    set(DEVICE "Mi100GPU" CACHE PATH "" FORCE)
  else()
    set(DEVICE "EPYCROME" CACHE PATH "" FORCE)
  endif()

endif()

if( ENABLE_CUDA )
  set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "")
  set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_ROOT} CACHE PATH "")
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "")
  set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
  set(extraComponentsLinkList ${extraComponentsLinkList} cuda )
  set(CMAKE_CUDA_STANDARD 20 CACHE STRING "")
  set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
endif()

include( ../blt/SetupBLT.cmake )

set (USE_VECTOR ON CACHE BOOL "" FORCE)

################################
# SEM USE OMP
################################
# the solver is using OpenMP+Vector on CPU
if(USE_OMP)
  set (SOLVER "OMP" CACHE PATH "" FORCE)
  message(STATUS "BUILDING SOLVER including OpenMP+Vector on CPU")
  set (USE_KOKKOS OFF CACHE BOOL "")
  set (USE_RAJA OFF CACHE BOOL "")
  set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (ENABLE_CUDA OFF CACHE STRING "" FORCE)
  set (USE_VECTOR ON CACHE BOOL "" FORCE)
  set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  message ("-- USE_OMP flag is set: USE_VECTOR= ${USE_VECTOR} ")
else()
    message ("-- USE_OMP flag is set to be false and OMP is not enabled")
endif()
################################


################################
# RAJA
################################
# the solver is using RAJA+LVARRAY+CUDA on GPU
if ( USE_RAJA)
   # disable OMP and KOKKOS
   set (USE_LVARRAY ON CACHE BOOL "" FORCE)
   set (USE_VECTOR OFF CACHE BOOL "" FORCE)
   set (USE_OMP OFF CACHE BOOL "")
   set (USE_KOKKOS OFF CACHE BOOL "")
   # set solver
   message( " " )
   message( STATUS "-------------------RAJA SETTINGS---------------------------------------" )
   set (SOLVER "Raja" CACHE PATH "" FORCE)
   message(STATUS "BUILDING SOLVER including RAJA ")

   set (RAJA_DIR "$ENV{RAJA_DIR}")
   if(DEFINED RAJA_DIR)
       message(STATUS "USE RAJA in RAJA_DIR = ${RAJA_DIR}")
       find_package(RAJA REQUIRED
                    PATHS ${RAJA_DIR}
                    NO_DEFAULT_PATH)
       get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(RAJA PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")
       message(STATUS "RAJA_INCLUDE_DIRS: ${RAJA_INCLUDE_DIRS}")
       set(ENABLE_RAJA ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} RAJA )
   else()
       message(FATAL_ERROR "ProxyApp requires RAJA, set RAJA_DIR to the RAJA installation directory.")
   endif()

   if( ENABLE_HIP )
     set(RAJA_ENABLE_HIP ON CACHE BOOL "" FORCE)
     set(RAJA_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
     include_directories("/opt/rocm-6.0.0/include")
   endif()

   ################################
   # Umpire
   ################################
   set (UMPIRE_DIR "$ENV{UMPIRE_DIR}")
   if(DEFINED UMPIRE_DIR)
       message(STATUS " USE RAJA: RAJA must be compiled with UMPIRE in UMPIRE_DIR = ${UMPIRE_DIR}")
       find_package(umpire REQUIRED
                    PATHS ${UMPIRE_DIR}
                    NO_DEFAULT_PATH)
       set(ENABLE_UMPIRE ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} umpire)
   else()
       message(FATAL_ERROR "ProxyApp requires Umpire, set UMPIRE_DIR to the Umpire installation directory.")
   endif()

   ################################
   # CHAI
   ################################
   set (CHAI_DIR "$ENV{CHAI_DIR}")
   if(DEFINED CHAI_DIR)
       message(STATUS " USE RAJA in CHAI_DIR = ${CHAI_DIR}")
       find_package(chai REQUIRED
                    PATHS ${CHAI_DIR}
                    NO_DEFAULT_PATH)
       get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(chai PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")
       set(ENABLE_CHAI ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} chai)
   else()
     message(FATAL_ERROR "ProxyApp requires CHAI, set CHAI_DIR to the CHAI installation directory.")
   endif()

   set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
   set (dependencyList lvarray RAJA chai)
   add_subdirectory( LvArray )
   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LvArray/src)

   if( ENABLE_HIP )
     set (dependencyList ${dependencyList} blt::hip blt::hip_runtime)
   endif()

   message( STATUS "-------------------END RAJA SETTINGS-----------------------------------" )
   message( " " )
else()
    message ("-- USE_RAJA flag is set to be false and RAJA is not enabled")
endif()
################################


################################
# KOKKOS
################################
# the solver is using KOKKOS+CUDA if USE_KOKKOS is ON
if (NOT USE_KOKKOS)
  message ("-- USE_KOKKOS flag is set to be false and KOKKOS is not enabled")
elseif (EXISTS $ENV{KOKKOS_DIR}/KokkosConfig.cmake)
  set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (USE_VECTOR OFF CACHE BOOL "" FORCE)
  message ("-- USE_KOKKOS flag is set to be true and KOKOS is now enabled")
  message("-- Found the KOKKOS library in $ENV{KOKKOS_DIR}")
  set (SOLVER "Kokkos" CACHE PATH "" FORCE)

  # KOKOS uses OMP
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

  # setup KOKOS
  set(Kokkos_ROOT "$ENV{KOKKOS_DIR}")
  find_package(Kokkos REQUIRED)
  set (Kokkos_INCLUDE_DIRS "${Kokkos_ROOT}/../../../include")
  include_directories(${Kokkos_INCLUDE_DIRS})

  set (kokkosList "-L${Kokkos_ROOT}/../.. -lkokkoscontainers -lkokkoscore -lkokkossimd ")
  set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  set (dependencyList ${dependencyList} ${kokkosList})
  message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
  if( ENABLE_HIP )
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
    set (dependencyList ${dependencyList} Kokkos::kokkos)
    message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
    include_directories("/opt/rocm-6.0.0/include")
  elseif ( ENABLE_CUDA )
    set (dependencyList "${dependencyList} -lcuda")
  endif()

elseif (NOT DEFINED $ENV{KOKKOS_DIR})
  set (USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not defined  not enabled")
else ()
  set (USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not the correct path and is not enabled")
endif()
################################

################################
## Shiva
###############################
#-------------------------------------------------------------------
# setup for Shiva
## -- the default is using Shiva library
## -- if USE_SHIVA=OFF it will use the SEM as in GEOS
#-------------------------------------------------------------------
if( USE_SHIVA )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Shiva/src/ShivaConfig.hpp.in
                  ${CMAKE_BINARY_DIR}/src/ShivaConfig.hpp )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Shiva/tpl/camp/include/camp/config.in.hpp
                  ${PROJECT_BINARY_DIR}/src/camp/config.hpp)
  message(STATUS "USE SHIVA as a submodule")
endif()
###############################

configure_file( utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/commonConfig.hpp)

add_subdirectory( sem )
add_subdirectory( fd )
#add_subdirectory( tests )
