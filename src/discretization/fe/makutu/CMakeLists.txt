set( finiteElement_headers
     include/LagrangeBasis1.h
     include/LagrangeBasis2.h
     include/LagrangeBasis3GL.h
     include/LagrangeBasis4GL.h
     include/LagrangeBasis5GL.h
     include/Qk_Hexahedron_Lagrange_GaussLobatto.h
     include/tensorops.h
)

add_library( finiteElement
             INTERFACE
             ${finiteElement_headers})

target_include_directories( finiteElement
                            INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
                            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                            $<INSTALL_INTERFACE:include> )

install( FILES ${finiteElement_headers}
         DESTINATION include/finiteElement )
