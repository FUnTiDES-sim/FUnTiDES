set( finiteElement_headers
     include/LagrangeBasis1.hpp
     include/LagrangeBasis2.hpp
     include/LagrangeBasis3GL.hpp
     include/LagrangeBasis4GL.hpp
     include/LagrangeBasis5GL.hpp
     include/Qk_Hexahedron_Lagrange_GaussLobatto.hpp
     include/tensorops.hpp
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
