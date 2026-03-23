add_library(funtides_inversion)

add_subdirectory(acoustic)
add_subdirectory(elastic)

set(INV_LIBS discretization proxy_utils proxy_model_struct proxy_model_unstruct funtides_inversion_api proxy_model_api)

target_link_kokkos_if_enabled(funtides_inversion)

target_link_libraries(funtides_inversion PUBLIC ${INV_LIBS})

install(TARGETS funtides_inversion
  EXPORT funtides_inversion_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
