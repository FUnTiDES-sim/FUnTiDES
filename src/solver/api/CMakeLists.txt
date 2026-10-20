add_library(fun_solver_api INTERFACE)
add_library(funtides::solver_api ALIAS fun_solver_api)
set_target_properties(fun_solver_api PROPERTIES EXPORT_NAME solver_api)

target_include_directories(fun_solver_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/funtides/solver>
)

target_link_libraries(fun_solver_api
  INTERFACE
    fun_core
    fun_model
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
