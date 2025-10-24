#-------------------------------------------------------------------
# Python wrapping configuration
#-------------------------------------------------------------------

if(NOT ENABLE_PYWRAP)
  message("-- ENABLE_PYWRAP is set to false. Python wrapping is not enabled")
elseif(ENABLE_PYWRAP)
  message(STATUS "-- Enabling Python wrapping libraries. Using pybind11")
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../external/pybind11)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE CACHE BOOL "Disable global LTO" FORCE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # Set Python install directory to avoid issues with conda and adios (default will try to install in /lib/pythonX.Y/site-packages)
  set(CMAKE_INSTALL_PYTHONDIR "${CMAKE_INSTALL_PREFIX}/python" CACHE PATH "Python install directory")
endif()