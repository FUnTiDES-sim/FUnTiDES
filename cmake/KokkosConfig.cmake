#-------------------------------------------------------------------
# Kokkos configuration
#-------------------------------------------------------------------

# the solver is using KOKKOS if USE_KOKKOS is ON
if(NOT USE_KOKKOS)
  message(STATUS "USE_KOKKOS flag is set to false. KOKKOS is not enabled")
else()
  message(STATUS "Building the application with Kokkos.")
  set(USE_VECTOR OFF CACHE BOOL "" FORCE)

  find_package(Kokkos REQUIRED)
endif()
