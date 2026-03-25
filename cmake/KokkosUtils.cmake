#-------------------------------------------------------------------
# Kokkos related functions
#-------------------------------------------------------------------

function(target_link_kokkos_if_enabled target_name)
    target_link_libraries(${target_name} PUBLIC Kokkos::kokkos)
    target_compile_definitions(${target_name} PUBLIC USE_KOKKOS)
endfunction()
