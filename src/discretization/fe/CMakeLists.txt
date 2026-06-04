add_library(fun_discretization INTERFACE)

target_include_directories( fun_discretization INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/makutu/include/>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/>
                          )
message(STATUS "Binary include dir: ${CMAKE_BINARY_DIR}/include")

target_link_libraries( fun_discretization INTERFACE fun_utils )
