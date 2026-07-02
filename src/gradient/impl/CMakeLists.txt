add_library(fun_gradient)


add_subdirectory(acoustic)
add_subdirectory(common)
add_subdirectory(elastic)

set(INV_LIBS fun_discretization fun_utils fun_model fun_gradient_api)

target_link_kokkos_if_enabled(fun_gradient)

target_link_libraries(fun_gradient PUBLIC ${INV_LIBS})

install(TARGETS fun_gradient
  EXPORT fun_gradient_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
