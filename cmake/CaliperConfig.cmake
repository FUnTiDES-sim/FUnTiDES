#-------------------------------------------------------------------
# Caliper configuration
#-------------------------------------------------------------------

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
  message("-- USE_CALIPER is enabled but unknown error occurred.")
  set(USE_CALIPER OFF CACHE BOOL "")
endif()