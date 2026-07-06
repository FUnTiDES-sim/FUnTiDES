add_library(fun_discretization_tensorial INTERFACE)

target_include_directories(fun_discretization_tensorial INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>)

target_link_libraries(fun_discretization_tensorial INTERFACE
                      fun_utils
                      fun_discretization_common)
