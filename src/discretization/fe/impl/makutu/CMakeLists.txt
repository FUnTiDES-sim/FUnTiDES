# makutu back-end (header only): Qk hexahedron kernels on Gauss-Lobatto-Legendre
# nodes. The 1D bases come from fun_discretization_common.
add_library(fun_discretization_makutu INTERFACE)

target_include_directories(fun_discretization_makutu INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
                          )

target_link_libraries(fun_discretization_makutu INTERFACE
                      fun_core
                      fun_discretization_common)
