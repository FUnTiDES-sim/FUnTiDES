add_library(fun_core INTERFACE)
add_library(funtides::core ALIAS fun_core)

target_include_directories(fun_core
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/core/include>
)
