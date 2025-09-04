add_library(proxy_model_struct
  src/model_struct.cc
)

target_include_directories(proxy_model_struct
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(proxy_model_struct
  PUBLIC
    proxy_model_api
)

target_link_kokkos_if_enabled(proxy_model_struct)