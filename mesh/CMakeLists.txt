add_library(mesh
    src/SEMmesh.cpp
)

target_include_directories(mesh PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

set(MESH_LIBS utils)

if (USE_KOKKOS)
    list(APPEND MESH_LIBS Kokkos::kokkos)
endif()

target_link_libraries(mesh PUBLIC ${MESH_LIBS})
