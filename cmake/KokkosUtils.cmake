#-------------------------------------------------------------------
# Kokkos related functions
#-------------------------------------------------------------------

function(target_link_kokkos_if_enabled target_name)
  if (USE_KOKKOS)
    target_link_libraries(${target_name} PUBLIC Kokkos::kokkos)
  endif()
endfunction()