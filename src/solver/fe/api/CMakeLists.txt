add_library(proxy_solver_api INTERFACE)

target_include_directories(proxy_solver_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(proxy_solver_api
  INTERFACE
    proxy_model_api
    proxy_utils
)