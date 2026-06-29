add_subdirectory(common)
add_subdirectory(makutu)

# Aggregate of every implementation target.
add_library(fun_discretization_impl INTERFACE)
target_link_libraries(fun_discretization_impl INTERFACE
                      fun_discretization_common
                      fun_discretization_makutu
                     )
