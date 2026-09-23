# Header-only code shared by the discretization back-ends: the 1D Lagrange
# bases, the math helpers, the type selectors of Integrals.h, and
# QkHexahedronBase (derived by no back-end yet, see docs/design-red-flags.md).
add_library(fun_discretization_common INTERFACE)

target_include_directories(fun_discretization_common INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
                          )

# common depends on the abstract layer (QkHexahedronBase derives from the tag,
# Integrals.h uses DiscretizationKind). It does NOT link makutu: the back-end
# headers it includes are resolved transitively through fun_discretization.
target_link_libraries(fun_discretization_common INTERFACE
                       fun_core
                       fun_discretization_api
                      )
