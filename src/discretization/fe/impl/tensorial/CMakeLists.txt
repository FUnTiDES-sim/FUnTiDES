add_library(fun_discretization_tensorial INTERFACE)
add_library(funtides::discretization_tensorial ALIAS fun_discretization_tensorial)
set_target_properties(fun_discretization_tensorial PROPERTIES
                      EXPORT_NAME discretization_tensorial)

target_include_directories(fun_discretization_tensorial INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
                           $<INSTALL_INTERFACE:include/funtides/discretization>)

target_link_libraries(fun_discretization_tensorial INTERFACE
                      fun_core
                      fun_discretization_common)

install_funtides_headers(discretization ${CMAKE_CURRENT_SOURCE_DIR}/include)
