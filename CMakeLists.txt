cmake_minimum_required( VERSION 3.22.1 )

#-----------------------------------------------------
# Cmake debug flag
#-----------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE OFF)

#--------------------------------------------------------
# Configure compiler/project
#--------------------------------------------------------
project(proxyapp LANGUAGES CXX)
include(cmake/CompilerConfig.cmake) # This must go *after* project()
include(cmake/ProxyConfig.cmake)    # This must go *before* other configurations

#--------------------------------------------------------
# Configure
#--------------------------------------------------------
include(cmake/PywrapConfig.cmake)   # This must go *before* add_subdirectory(kokkos)
include(cmake/KokkosConfig.cmake)
include(cmake/ADIOS2.cmake)
include(cmake/Googletest.cmake)

#--------------------------------------------------------
# Custom cmake functions
#--------------------------------------------------------
include(cmake/KokkosUtils.cmake)
include(cmake/ProxyUtils.cmake)
include(cmake/PywrapUtils.cmake)

#--------------------------------------------------------
# Subdirectories
#--------------------------------------------------------
add_subdirectory(src)
add_subdirectory(tests)

#--------------------------------------------------------
# Print configuration summary
#--------------------------------------------------------
print_configuration_summary()
