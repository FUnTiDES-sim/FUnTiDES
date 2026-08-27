set(MODEL_TEST_SOURCES
  impl/builder/cartesian/test_cartesian_struct_boundary_classifier.cc
  impl/builder/cartesian/test_cartesian_struct_builder.cc
  impl/builder/cartesian/test_cartesian_unstruct_boundary_classifier.cc
  impl/builder/cartesian/test_cartesian_model_file_reader.cc
  impl/builder/cartesian/test_cartesian_unstruct_builder.cc
  impl/mesh/test_topology_factory.cc
  impl/mesh/faceconnectivity/test_face_connectivity_unstruct.cpp
  impl/mesh/faceconnectivity/test_face_connectivity_struct.cpp
  impl/mesh/struct/test_struct_mesh.cc
  impl/mesh/struct/test_gllpoints.cc
  impl/mesh/unstruct/test_unstruct_mesh.cc
)

add_gtest_group(model_tests
  SOURCES ${MODEL_TEST_SOURCES}
  EXTRA_LINKS fun_model
)

# Kept as its own executable: defines the same TopologyFactoryTest fixture/tests
# as impl/mesh/test_topology_factory.cc (duplicate content) — would collide with
# it if merged into model_tests.
add_gtest(cartesian_partitioner_test
  impl/builder/cartesian/test_cartesian_partitioner.cc
  fun_model
)
