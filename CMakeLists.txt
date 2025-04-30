cmake_minimum_required( VERSION 3.22.1 )


#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------
set (CMAKE_VERBOSE_MAKEFILE ON )

option (COMPILE_SEM "Compile Spectral Elements Method simulation" ON)
option (COMPILE_FD "Compile finit elements simulation" ON)
#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------
# use of Caliper to instrument code
option (USE_CALIPER "Use Caliper instrumentation" OFF)

# use of EZV to instrument code
option (USE_EZV "Use EZV instrumentation" OFF)

# use OMP to parallelise loops
option (USE_OMP "Use OpenMP to parallelise loops" OFF)

# use KOKOS to expose parallelism loop
option (USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)

# Enable or not cuda
option (ENABLE_CUDA "Enable cuda compilation" OFF)

# Enable the python wrapping compilation
option (ENABLE_PYWRAP "Enable python binding compilation with pybind11" OFF)

# Chose integral version
option (USE_SEMCLASSIC "Use CLASSIC " OFF)
option (USE_SEMOPTIM "Use OPTIM " ON)
option (USE_SHIVA "Use SHIVA " OFF)

# Chose data model
option (USE_SIMPLE_MESH "Use Simple mesh" ON)

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

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
# So make install will copy pykokkos-base onto proxy folder
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "." CACHE PATH "Install path prefix" FORCE)
endif()
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

 # CC options
set(CMAKE_C_COMPILER  /usr/bin/gcc-12  CACHE STRING "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
# C++ options
set(CMAKE_CXX_COMPILER /usr/bin/g++-12  CACHE PATH "")
message(STATUS "--CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG " CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable gnu++ extensions" FORCE)

project( proxyapp )

################################
# SEM USE OMP
################################
# the solver is using OpenMP+Vector on CPU
if(USE_OMP AND NOT USE_KOKKOS)
  message(STATUS "BUILDING SOLVER including OpenMP+Vector on CPU")
  find_package(OpenMP REQUIRED)
  set (ENABLE_CUDA OFF CACHE STRING "" FORCE)
  set (USE_VECTOR ON CACHE BOOL "" FORCE)
  message ("-- USE_OMP flag is set: USE_VECTOR= ${USE_VECTOR} ")
else()
    message ("-- USE_OMP flag is set to be false or Kokkos activated. OMP is not enabled")
endif()

################################
# PYTHON WRAPPING
################################
if (ENABLE_PYWRAP)
  message(STATUS "-- Enabling Python wrapping librairies. Using pybind11")
  find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE CACHE BOOL "Disable global LTO" FORCE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

################################
# KOKKOS
################################
# the solver is using KOKKOS+CUDA if USE_KOKKOS is ON
if (NOT USE_KOKKOS)
  message ("-- USE_KOKKOS flag is set to false. KOKKOS is not enabled")
else ()
  message (STATUS "Building the application with Kokkos.")
  set(Kokkos_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/kokkos)
  set(Kokkos_BINARY_DIR ${CMAKE_BINARY_DIR}/external/kokkos/)

  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_SHARED ON CACHE BOOL "" FORCE)
  if (ENABLE_CUDA)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_UVM ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
    message (STATUS "-- Activating CUDA with Kokkos.")
  endif()
  add_subdirectory(${Kokkos_SOURCE_DIR})
  message (STATUS "-- Kokkos Activated. Source are at ${Kokkos_SOURCE_DIR}.")
  if(ENABLE_PYWRAP)
    set(ENABLE_INTERNAL_KOKKOS OFF CACHE BOOL "" FORCE)
    set(Kokkos_DIR "${Kokkos_BINARY_DIR}" CACHE PATH "Use external Kokkos build")
    message(STATUS "-- Enabling pykokkos base.")
    add_subdirectory(external/pykokkos-base)
    # Installing pykokkos-base python lib.
    set(PYKOKKOS_BUILD_DIR "${CMAKE_BINARY_DIR}/external/pykokkos-base/kokkos")
    file(GLOB PYKOKKOS_MODULE_SO "${PYKOKKOS_BUILD_DIR}/libpykokkos*.so")

      install(
        FILES ${PYKOKKOS_MODULE_SO}
        DESTINATION proxy
      )
  endif()
endif()

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

   set (dependencyList ${dependencyList} sdl2)
elseif(NOT DEFINED $ENV{EASYPAP_DIR})
    message("-- USE_EZV is enabled but EASYPAP_DIR is not defined.")
    set(USE_EZV OFF CACHE BOOL "")
else()
    message("-- USE_EZV flag is enabled but unknown error occured.")
    set(USE_EZV OFF CACHE BOOL "")
endif()


configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/include/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/utils/include/commonConfig.hpp)

 add_subdirectory(src)
