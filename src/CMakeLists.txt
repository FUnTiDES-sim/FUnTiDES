#-------------------------------------------------------------------
# The Solver Implementation is with Programming Models:
#
#  Users can choose from "OMP" or "KOKKOS" as following:
#
#  (0) The default is using sequential mode
#  (1) If USE_OMP is ON, it is using OpenMP+Vector on CPU.
#  (2) If USE_KOKKOS is ON, it is using KOKKOS+CUDA.
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
  set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

  # C++ options
  set(CMAKE_CXX_COMPILER /usr/bin/g++-12  CACHE PATH "")
  message( STATUS, "-----------------------CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
  set (CMAKE_CXX_FLAGS "-O3 -std=c++20 -fconcepts " CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG " CACHE STRING "")
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
    set(CMAKE_CUDA_FLAGS "-restrict --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -std=c++20 -fconcepts" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELEASE "-g -Xcompiler -DNDEBUG" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -Xcompiler -O0" CACHE STRING "")
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
  # set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (ENABLE_CUDA OFF CACHE STRING "" FORCE)
  set (USE_VECTOR ON CACHE BOOL "" FORCE)
  set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  message ("-- USE_OMP flag is set: USE_VECTOR= ${USE_VECTOR} ")
else()
    message ("-- USE_OMP flag is set to be false and OMP is not enabled")
endif()
################################

################################
# KOKKOS
################################
# the solver is using KOKKOS+CUDA if USE_KOKKOS is ON
if (NOT USE_KOKKOS)
  message ("-- USE_KOKKOS flag is set to be false and KOKKOS is not enabled")
elseif (EXISTS $ENV{KOKKOS_DIR}/KokkosConfig.cmake)
  # set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
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


################################
## Caliper
###############################
if(NOT USE_CALIPER)
    message("-- USE_CALIPER is set to false. CALIPER is not enabled")
elseif(EXISTS $ENV{CALIPER_DIR}/share/cmake/caliper/caliper-config.cmake)
    message("-- Using Caliper for code instrumentation")
    # Manually set the paths
    set(CALIPER_INCLUDE_DIRS "$ENV{CALIPER_DIR}/include")
    set(CALIPER_LIBRARIES "$ENV{CALIPER_DIR}/lib/libcaliper.so")
    message(STATUS "Caliper Include Dirs: ${CALIPER_INCLUDE_DIRS}")
    message(STATUS "Caliper Libraries: ${CALIPER_LIBRARIES}")
    include_directories(${CALIPER_INCLUDE_DIRS})
    link_directories(${CALIPER_DIR}/lib)
    # Find package
    set(caliper_DIR "$ENV{CALIPER_DIR}/share/cmake/caliper")
    find_package(caliper REQUIRED)
    message("-- Caliper founded at $ENV{CALIPER_DIR}")
    
    # Register Caliper with BLT (if using BLT)
    blt_register_library(
	NAME caliper
	INCLUDES ${CALIPER_INCLUDE_DIRS}
	LIBRARIES ${CALIPER_LIBRARIES}
    )
    set (dependencyList ${dependencyList} caliper)
        
    add_definitions(-DUSE_CALIPER)
elseif(NOT DEFINED $ENV{CALIPER_DIR})
    message("-- USE_CALIPER is enabled but CALIPER_DIR is not defined.")
    set(USE_CALIPER OFF CACHE BOOL "")
else()
    message("-- USE_CALIPER is enabled but unknown error occured.")
    set(USE_CALIPER OFF CACHE BOOL "")
endif()

################################
## Easy Viz
###############################
if(NOT USE_EZV)
    message("-- USE_EZV is set to false. EZV is not enabled")
elseif(EXISTS $ENV{EASYPAP_DIR}/lib/cmake/Easypap/ezvConfig.cmake)
    message("-- Using EZV for code instrumentation")
    # Find Scotch dependancy
    message("-- Scotch is required for ezv to work.")
    message("-- SCOTCH_DIR variable must be set at scotch root directory.")
    set(SCOTCH_DIR "$ENV{SCOTCH_DIR}/lib/cmake/scotch")
    message("-- SCOTCH_DIR now set at $ENV{SCOTCH_DIR}.")
    # Find package
    set(Easypap_ROOT "$ENV{EASYPAP_DIR}/lib/cmake/Easypap")
    find_package(Easypap COMPONENTS ezv REQUIRED)
    set(EZV_INCLUDE_DIRS "$ENV{EASYPAP_DIR}/include")
    set(EZV_LIBRARIES "$ENV{EASYPAP_DIR}/lib/libezv.so")
    message("-- EZV founded at $ENV{EASYPAP_DIR}")
    message(STATUS "EZV Include Dirs: ${EZV_INCLUDE_DIRS}")
    message(STATUS "EZV Libraries: ${EZV_LIBRARIES}")
    # SDL 2 configuration
    message("-- SLD2 must be installed for EZV.")
    find_package(SDL2 REQUIRED)
    include_directories(${SDL2_INCLUDE_DIRS})

    # Register EZM with BLT (if using BLT)
    blt_register_library(
	NAME ezv
	INCLUDES ${EZV_INCLUDE_DIRS}
	LIBRARIES ${EZV_LIBRARIES}
    )
    set (dependencyList ${dependencyList} ezv)
    add_definitions(-DUSE_EZV)
    # Register SDL2 with BLT (if using BLT)
    blt_register_library(
	NAME sdl2
	INCLUDES ${SDL2_INCLUDE_DIRS}
	LIBRARIES ${SDL2_LIBRARIES}
    )
    set (dependencyList ${dependencyList} sdl2)
elseif(NOT DEFINED $ENV{EASYPAP_DIR})
    message("-- USE_EZV is enabled but EASYPAP_DIR is not defined.")
    set(USE_EZV OFF CACHE BOOL "")
else()
    message("-- USE_EZV is enabled but unknown error occured.")
    set(USE_EZV OFF CACHE BOOL "")
endif()



configure_file( utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/commonConfig.hpp)

add_subdirectory( sem )
add_subdirectory( fd )
#add_subdirectory( tests )
