add_library(cartesian_mesh
    src/cartesian_sem_mesh.cc
)

target_include_directories(cartesian_mesh PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../meshAPI/include>
)

set(MESH_LIBS utils)

if (USE_KOKKOS)
    list(APPEND MESH_LIBS Kokkos::kokkos)
endif()

target_link_libraries(cartesian_mesh PUBLIC ${MESH_LIBS})


# Adding test case
add_executable(test_cartesian src/test_mesh.cc)
target_include_directories(test_cartesian PRIVATE  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_cartesian
    PRIVATE
      cartesian_mesh
      utils
  )
