add_library(funtides_gradient_api INTERFACE)

target_include_directories(funtides_gradient_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(funtides_gradient_api
  INTERFACE
    proxy_model_api
    fun_utils
)
