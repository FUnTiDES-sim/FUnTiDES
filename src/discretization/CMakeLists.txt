add_library(discretization INTERFACE)

target_include_directories(discretization
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe/SEMKernels/src/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe/SEMKernels/src/finiteElement/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe>
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fd>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    ${CMAKE_BINARY_DIR}/include
)
message(STATUS "Binary include dir: ${CMAKE_BINARY_DIR}/include")

add_subdirectory( fe/SEMKernels )

target_link_libraries( discretization INTERFACE proxy_utils )
# target_link_libraries( discretization INTERFACE Shiva::shiva )
