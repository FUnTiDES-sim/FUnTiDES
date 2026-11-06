cmake_minimum_required( VERSION 3.22.1 )

#-----------------------------------------------------
# Cmake debug flag
#-----------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE OFF)

#--------------------------------------------------------
# Configure compiler/project
#--------------------------------------------------------
project(funtides VERSION 0.2.0 LANGUAGES CXX)
include(cmake/CompilerConfig.cmake) # This must go *after* project()
include(cmake/ProxyConfig.cmake)    # This must go *before* other configurations

#--------------------------------------------------------
# Configure
#--------------------------------------------------------
include(cmake/GenerateVariables.cmake)
include(cmake/PywrapConfig.cmake)   # This must go *before* kokkos configuration
include(cmake/MPIConfig.cmake)
if(ENABLE_TESTS)
  include(cmake/GoogleTestConfig.cmake)
  include(cmake/GoogleBenchConfig.cmake)
endif()

#--------------------------------------------------------
# Custom cmake functions
#--------------------------------------------------------
include(cmake/KokkosUtils.cmake)
include(cmake/ProxyUtils.cmake)
include(cmake/PywrapUtils.cmake)
if(ENABLE_TESTS)
  include(cmake/GoogleBenchUtils.cmake)
  include(cmake/GoogleTestUtils.cmake)
endif()

#--------------------------------------------------------
# Subdirectories
#--------------------------------------------------------
add_subdirectory(src)
if(ENABLE_TESTS)
  add_subdirectory(tests)
endif()

#--------------------------------------------------------
# Print configuration summary
#--------------------------------------------------------
print_configuration_summary()
