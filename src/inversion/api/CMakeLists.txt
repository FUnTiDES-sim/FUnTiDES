add_library(funtides_inversion_api INTERFACE)

target_include_directories(funtides_inversion_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(funtides_inversion_api
  INTERFACE
    proxy_model_api
    proxy_utils
)
