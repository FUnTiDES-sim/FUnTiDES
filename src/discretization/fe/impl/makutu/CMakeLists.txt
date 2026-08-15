# Concrete back-ends (Lagrange-GaussLobatto + tensorial GEMM) and the GLL bases.
add_library(fun_discretization_makutu INTERFACE)

target_include_directories(fun_discretization_makutu INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
                          )

target_link_libraries(fun_discretization_makutu INTERFACE
                      fun_utils
                      fun_discretization_common)
