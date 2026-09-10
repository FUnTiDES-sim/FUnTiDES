#-------------------------------------------------------------------
# Proxy app configuration
#-------------------------------------------------------------------

# Discretization
option(COMPILE_SEM "Compile Spectral Elements Method simulation" ON)
option(COMPILE_DG "Compile Discontinuous Galerkin simulation" ON)
option(COMPILE_DG_SEM "Compile Discontinuous Galerkin coupled with Spectral Elements Method simulation" ON)
option(COMPILE_DG_PADAPTIVE "Compile p-adaptive Discontinuous Galerkin simulation" ON)

if (COMPILE_DG_SEM)
  set(COMPILE_DG ON)
  set(COMPILE_SEM ON)
  message(STATUS "DGSEM solver will be compiled, so DG and SEM solvers will also be compiled")
endif()

if(COMPILE_DG_PADAPTIVE)
  set(COMPILE_DG ON)
  message(STATUS "p-adaptive DG solver will be compiled, so DG solver will also be compiled")
endif()

# Programming models
option(USE_MPI "Enable MPI compilation" OFF)

# Python wrapping
option(ENABLE_PYWRAP "Enable python binding compilation with pybind11" OFF)

# Debugging options
option(PRINT_ALLOC_INFO "Printout memory allocation info" OFF)
# Build options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Install options
# So make install will copy pykokkos onto proxy folder
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "." CACHE PATH "Install path prefix" FORCE)
endif()

# Macro definitions
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/utils/include/common_config.h.in
               ${CMAKE_BINARY_DIR}/src/utils/include/common_config.h)
