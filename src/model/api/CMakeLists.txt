add_library(proxy_model_api INTERFACE)

target_include_directories(proxy_model_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(proxy_model_api
  INTERFACE
    proxy_utils
)