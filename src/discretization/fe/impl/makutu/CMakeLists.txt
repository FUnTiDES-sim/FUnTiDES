# makutu back-end (header only): Qk hexahedron kernels on Gauss-Lobatto-Legendre
# nodes. The 1D bases come from fun_discretization_common.
add_library(fun_discretization_makutu INTERFACE)
add_library(funtides::discretization_makutu ALIAS fun_discretization_makutu)
set_target_properties(fun_discretization_makutu PROPERTIES
                      EXPORT_NAME discretization_makutu)

target_include_directories(fun_discretization_makutu INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
                           $<INSTALL_INTERFACE:include/funtides/discretization>
                          )

target_link_libraries(fun_discretization_makutu INTERFACE
                      fun_core
                      fun_discretization_common)

install_funtides_headers(discretization ${CMAKE_CURRENT_SOURCE_DIR}/include)
