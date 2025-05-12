add_library(discretization INTERFACE)

target_include_directories(discretization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Shiva/src/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fd>
)
