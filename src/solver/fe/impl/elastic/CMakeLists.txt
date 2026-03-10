add_subdirectory(src/makutu)

target_include_directories(proxy_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
