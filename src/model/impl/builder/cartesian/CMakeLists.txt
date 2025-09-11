add_library(proxy_model_builder_cartesian
  src/cartesian_struct_builder.cc
  src/cartesian_unstruct_builder.cc
)

target_include_directories(proxy_model_builder_cartesian PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(proxy_model_builder_cartesian PUBLIC
  proxy_model_struct
  proxy_model_unstruct
)

target_link_kokkos_if_enabled(proxy_model_builder_cartesian)

install(TARGETS proxy_model_builder_cartesian
  EXPORT proxy_model_builder_cartesian_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)