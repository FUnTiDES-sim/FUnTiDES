add_library(discretization INTERFACE)

target_include_directories( discretization INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/makutu/include/>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/>
                          )
message(STATUS "Binary include dir: ${CMAKE_BINARY_DIR}/include")

target_link_libraries( discretization INTERFACE proxy_utils )
