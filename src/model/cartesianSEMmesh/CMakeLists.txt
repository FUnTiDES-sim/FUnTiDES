add_library(mesh
    src/cartesianSEMmesh.cpp
)

target_include_directories(mesh PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../meshAPI/include>
)

set(MESH_LIBS utils)

if (USE_KOKKOS)
    list(APPEND MESH_LIBS Kokkos::kokkos)
endif()

target_link_libraries(mesh PUBLIC ${MESH_LIBS})


# Adding test case
add_executable(testmesh src/test_mesh.cpp)
target_include_directories(testmesh PRIVATE  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(testmesh
    PRIVATE
      mesh
      utils
  )
