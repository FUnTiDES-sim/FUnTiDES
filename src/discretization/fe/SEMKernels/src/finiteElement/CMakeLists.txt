set( finiteElement_headers
     makutu/LagrangeBasis1.hpp
     makutu/LagrangeBasis2.hpp
     makutu/LagrangeBasis3GL.hpp
     makutu/LagrangeBasis4GL.hpp
     makutu/LagrangeBasis5GL.hpp
     makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp
     makutu/tensorops.hpp
)

if( ENABLE_Shiva )
    list( APPEND finiteElement_headers
          shiva/SEMQkGLIntegralsShiva.hpp
    )
endif()

set( finiteElement_dependencies
     shiva
     )

add_library( finiteElement
             INTERFACE
             ${finiteElement_headers}
             ${finiteElement_dependencies} )

target_include_directories( finiteElement
                            INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
                            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                            $<INSTALL_INTERFACE:include> )

install( FILES ${finiteElement_headers}
         DESTINATION include/finiteElement )

if( SEMKERNELS_ENABLE_UNIT_TESTS )
    add_subdirectory( unitTests )
endif()
