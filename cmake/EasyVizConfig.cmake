#-------------------------------------------------------------------
# EasyViz configuration
#-------------------------------------------------------------------

if(NOT USE_EZV)
  message("-- USE_EZV is set to false. EZV is not enabled")
elseif(EXISTS $ENV{EASYPAP_DIR}/lib/cmake/Easypap/ezvConfig.cmake)
  message("-- Using EZV for code instrumentation")
  # Find Scotch dependency
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
  message("-- SDL2 must be installed for EZV.")
  find_package(SDL2 REQUIRED)
  include_directories(${SDL2_INCLUDE_DIRS})

  set(dependencyList ${dependencyList} sdl2)
elseif(NOT DEFINED $ENV{EASYPAP_DIR})
  message("-- USE_EZV is enabled but EASYPAP_DIR is not defined.")
  set(USE_EZV OFF CACHE BOOL "")
else()
  message("-- USE_EZV flag is enabled but unknown error occurred.")
  set(USE_EZV OFF CACHE BOOL "")
endif()