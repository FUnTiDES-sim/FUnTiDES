# Header-only interface target: discretization contract and back-end selector.
add_library(fun_discretization_api INTERFACE)

target_include_directories(fun_discretization_api INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
                          )
