# Header-only interface target: discretization contract and back-end selector.
add_library(fun_discretization_api INTERFACE)
add_library(funtides::discretization_api ALIAS fun_discretization_api)
set_target_properties(fun_discretization_api PROPERTIES
                      EXPORT_NAME discretization_api)

target_include_directories(fun_discretization_api INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
                           $<INSTALL_INTERFACE:include/funtides/discretization>
                          )

install_funtides_headers(discretization ${CMAKE_CURRENT_SOURCE_DIR})
