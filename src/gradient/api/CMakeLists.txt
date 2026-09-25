add_library(fun_gradient_api INTERFACE)
add_library(funtides::gradient_api ALIAS fun_gradient_api)
set_target_properties(fun_gradient_api PROPERTIES EXPORT_NAME gradient_api)

target_include_directories(fun_gradient_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/funtides/gradient>
)

target_link_libraries(fun_gradient_api
  INTERFACE
    fun_model
    fun_core
)

install_funtides_headers(gradient ${CMAKE_CURRENT_SOURCE_DIR}/include)
