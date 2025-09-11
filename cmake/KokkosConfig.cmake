#-------------------------------------------------------------------
# Kokkos configuration
#-------------------------------------------------------------------

# the solver is using KOKKOS if USE_KOKKOS is ON
if(NOT USE_KOKKOS)
  message(STATUS "USE_KOKKOS flag is set to false. KOKKOS is not enabled")
else()
  message(STATUS "Building the application with Kokkos.")
  set(USE_VECTOR OFF CACHE BOOL "" FORCE)

  set(Kokkos_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../external/kokkos)
  set(Kokkos_BINARY_DIR ${CMAKE_BINARY_DIR}/external/kokkos/)

  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_SHARED ON CACHE BOOL "" FORCE)

  if(ENABLE_CUDA)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_UVM ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_RELOCATE_DEVICE_CODE OFF CACHE BOOL "" FORCE)
    message(STATUS "Activating CUDA with Kokkos.")
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(STATUS "Activating Kokkos debug mode")
    set(Kokkos_ENABLE_DEBUG ON CACHE BOOL "" FORCE)
  endif()

  add_subdirectory(${Kokkos_SOURCE_DIR})
  message(STATUS "Kokkos Activated. Source are at ${Kokkos_SOURCE_DIR}.")

  if(ENABLE_PYWRAP)
    # Enabling pykokkos-base and reused Kokkos build above
    set(ENABLE_MEMORY_TRAITS OFF)
    set(Kokkos_ENABLE_PYTHON ON)
    set(Kokkos_ENABLE_DEBUG OFF)
    set(BUILD_SHARE_LIBS ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
    set(ENABLE_INTERNAL_KOKKOS OFF CACHE BOOL "" FORCE)
    set(ENABLE_INTERNAL_PYBIND11 OFF CACHE BOOL "" FORCE)
    set(Kokkos_DIR "${Kokkos_BINARY_DIR}" CACHE PATH "Use external Kokkos build")
    message(STATUS "Enabling pykokkos base.")
    add_subdirectory(external/pykokkos-base)

    # Installing pykokkos-base python lib
    set(PYKOKKOS_BUILD_DIR "${CMAKE_BINARY_DIR}/external/pykokkos-base/kokkos")
    file(GLOB PYKOKKOS_MODULE_SO "${PYKOKKOS_BUILD_DIR}/libpykokkos*.so")
    message(STATUS "Installing pykokkos module: ${PYKOKKOS_MODULE_SO}")
    install(
      FILES ${PYKOKKOS_MODULE_SO}
      DESTINATION lib64
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
  endif()
endif()