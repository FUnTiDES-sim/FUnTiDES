#-------------------------------------------------------------------
# Proxy app configuration
#-------------------------------------------------------------------

# Discretization
option(COMPILE_SEM "Compile Spectral Elements Method simulation" ON)
option(COMPILE_FD "Compile finite elements simulation" ON)

# Programming models
option(USE_VECTOR "Use vectors." OFF)
option(USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)
option(USE_KOKKOS_TEAMS "use hierarchical parallelism in Kokkos" OFF)
option(ENABLE_CUDA "Enable cuda compilation" OFF)

# Python wrapping
option(ENABLE_PYWRAP "Enable python binding compilation with pybind11" OFF)

# Debugging options
option(FD_SAVE_SNAPSHOTS "Save snapshots for FD-proxy" OFF)
option(PRINT_ALLOC_INFO "Printout memory allocation info" OFF)

# Build options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Install options
# So make install will copy pykokkos-base onto proxy folder
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "." CACHE PATH "Install path prefix" FORCE)
endif()

# Macro definitions
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/utils/include/common_config.h.in
               ${CMAKE_BINARY_DIR}/src/utils/include/common_config.h)
