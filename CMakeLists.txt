cmake_minimum_required( VERSION 3.22.1 )

#-----------------------------------------------------
# Cmake debug flag
#-----------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE OFF)

#--------------------------------------------------------
# Configure compiler/project
#--------------------------------------------------------
include(cmake/CompilerConfig.cmake) # This must go *before* project()
project(proxyapp)
include(cmake/ProxyConfig.cmake)    # This must go *before* other configurations

#--------------------------------------------------------
# Configure
#--------------------------------------------------------
include(cmake/OpenmpConfig.cmake)
include(cmake/PywrapConfig.cmake)   # This must go *before* add_subdirectory(kokkos)
include(cmake/KokkosConfig.cmake)
include(cmake/CaliperConfig.cmake)
include(cmake/EasyVizConfig.cmake)

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

#--------------------------------------------------------
# Print configuration summary
#--------------------------------------------------------
print_configuration_summary()