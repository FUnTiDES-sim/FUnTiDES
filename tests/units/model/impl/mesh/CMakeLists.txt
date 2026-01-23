add_subdirectory(struct)
add_subdirectory(unstruct)

add_gtest(topology_factory_test
  test_topology_factory.cc
  proxy_model_api
)
