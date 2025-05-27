add_library(discretization INTERFACE)

target_include_directories( discretization INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe/SEMKernels/src/>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe/SEMKernels/src/finiteElement/>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe/SEMKernels/src/Shiva/src/>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fe>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fd> 
                            ${CMAKE_BINARY_DIR}/include
                            )
message(STATUS "Binary include dir: ${CMAKE_BINARY_DIR}/include")

include( fe/SEMKernels/cmake/blt/SetupBLT.cmake )
add_subdirectory( fe/SEMKernels )

add_dependencies( discretization finiteElements shiva)

