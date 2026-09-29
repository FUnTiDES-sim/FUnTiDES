add_subdirectory(common)
add_subdirectory(makutu)
add_subdirectory(tensorial)

# Aggregate of every implementation target.
add_library(fun_discretization_impl INTERFACE)
add_library(funtides::discretization_impl ALIAS fun_discretization_impl)
set_target_properties(fun_discretization_impl PROPERTIES
                      EXPORT_NAME discretization_impl)
target_link_libraries(fun_discretization_impl INTERFACE
                      fun_discretization_common
                      fun_discretization_makutu
                      fun_discretization_tensorial
                     )
