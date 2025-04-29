add_library(discretization INTERFACE)

target_include_directories(discretization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/discretization/include/discretization/fe>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/discretization/include/discretization/fd>
)
