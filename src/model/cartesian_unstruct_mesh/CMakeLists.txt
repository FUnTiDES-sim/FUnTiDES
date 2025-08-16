add_library(cartesianUnstructMesh
    src/cartesian_unstruct_mesh.cc
)

target_include_directories(cartesianUnstructMesh PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../meshAPI/include>
)

set(MESH_LIBS utils)

if (USE_KOKKOS)
    list(APPEND MESH_LIBS Kokkos::kokkos)
endif()

target_link_libraries(cartesianUnstructMesh PUBLIC ${MESH_LIBS})


# Adding test case
add_executable(testUmesh src/test_mesh.cc)
target_include_directories(testUmesh PRIVATE  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(testUmesh
    PRIVATE
      cartesianUnstructMesh
      utils
  )
