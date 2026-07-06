# Shared implementation: mother class QkHexahedronBase, DiscretizationTraits,
# and the compile-time / math helpers shared by every back-end.
add_library(fun_discretization_common INTERFACE)

target_include_directories(fun_discretization_common INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
                          )

# common depends on the abstract layer (QkHexahedronBase derives from the tag,
# Integrals.h uses DiscretizationKind). It does NOT link makutu: the back-end
# headers it includes are resolved transitively through fun_discretization.
target_link_libraries(fun_discretization_common INTERFACE
                       fun_utils
                       fun_discretization_api
                      )
