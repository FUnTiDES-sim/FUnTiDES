add_library(funtides_gradient)


add_subdirectory(acoustic)
add_subdirectory(common)
add_subdirectory(elastic)

set(INV_LIBS fun_discretization fun_utils proxy_model_struct proxy_model_unstruct funtides_gradient_api proxy_model_api)

target_link_kokkos_if_enabled(funtides_gradient)

target_link_libraries(funtides_gradient PUBLIC ${INV_LIBS})

install(TARGETS funtides_gradient
  EXPORT funtides_gradient_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
