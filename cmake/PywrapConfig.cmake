#-------------------------------------------------------------------
# Python wrapping configuration
#-------------------------------------------------------------------
# fix warning yield when install path is ./python
if(POLICY CMP0177)
    cmake_policy(SET CMP0177 NEW)
endif()

if(NOT ENABLE_PYWRAP)
  message("-- ENABLE_PYWRAP is set to false. Python wrapping is not enabled")
elseif(ENABLE_PYWRAP)
  message(STATUS "-- Enabling Python wrapping libraries. Using pybind11")
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  find_package(pybind11 REQUIRED)
  # add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../external/pybind11)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE CACHE BOOL "Disable global LTO" FORCE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # Set Python install directory to avoid issues with conda and adios (default will try to install in /lib/pythonX.Y/site-packages)
  # Relative to CMAKE_INSTALL_PREFIX, so it follows a later change of the prefix
  set(CMAKE_INSTALL_PYTHONDIR "python" CACHE STRING
      "Python install directory (relative to CMAKE_INSTALL_PREFIX, or absolute)")
endif()
