add_subdirectory(api)
add_subdirectory(impl)

# Umbrella target kept for backward compatibility: the rest of the project links
# `fun_discretization` and transitively gets the api + impl include paths.
add_library(fun_discretization INTERFACE)

target_link_libraries(fun_discretization INTERFACE
                      fun_utils
                      fun_discretization_api
                      fun_discretization_impl
                     )

message(STATUS "fun_discretization: api + impl{common,makutu} configured")
