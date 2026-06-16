add_subdirectory(src/makutu)

target_include_directories(fun_gradient
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
