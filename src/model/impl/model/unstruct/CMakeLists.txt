add_library(proxy_model_unstruct
  src/model_unstruct.cc
)

target_include_directories(proxy_model_unstruct
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(proxy_model_unstruct
  PUBLIC
    proxy_model_api
)

target_link_kokkos_if_enabled(proxy_model_unstruct)