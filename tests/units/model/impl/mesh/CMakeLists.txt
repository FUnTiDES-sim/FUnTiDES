add_subdirectory(struct)
add_subdirectory(unstruct)
add_subdirectory(faceconnectivity)

add_gtest(topology_factory_test
  test_topology_factory.cc
  fun_model
)
