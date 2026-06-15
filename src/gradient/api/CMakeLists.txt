add_library(fun_gradient_api INTERFACE)

target_include_directories(fun_gradient_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(fun_gradient_api
  INTERFACE
    fun_model
    fun_utils
)
