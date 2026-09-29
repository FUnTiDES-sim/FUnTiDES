target_include_directories(fun_solver_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
