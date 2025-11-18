add_library(proxy_model_builder_sep
  src/sep_builder.cc
)

target_include_directories(proxy_model_builder_sep PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(proxy_model_builder_sep PUBLIC
  proxy_model_struct
  proxy_model_unstruct
)

target_link_kokkos_if_enabled(proxy_model_builder_sep)

install(TARGETS proxy_model_builder_sep
  EXPORT proxy_model_builder_sep_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
