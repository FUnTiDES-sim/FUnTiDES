set(MODEL_TEST_SOURCES
  impl/builder/cartesian/test_cartesian_struct_boundary_classifier.cc
  impl/builder/cartesian/test_cartesian_struct_builder.cc
  impl/builder/cartesian/test_cartesian_unstruct_boundary_classifier.cc
  impl/builder/cartesian/test_cartesian_model_file_reader.cc
  impl/builder/cartesian/test_cartesian_unstruct_builder.cc
  impl/common/test_elasticity_utils.cc
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
# as tests/units/parallel/test_topology_factory.cc (duplicate content). The
# collision no longer applies since that file moved to parallel_tests, so this
# target could now be folded into model_tests — or the duplication removed.
add_gtest(cartesian_partitioner_test
  impl/builder/cartesian/test_cartesian_partitioner.cc
  "fun_model;fun_parallel"
)

# 2-D Cartesian builder and ModelUnstruct of dimension 2 (COMPILE_SEM_2D)
if(COMPILE_SEM_2D)
  add_gtest(gmsh_2d_reader_test
    impl/builder/gmsh/test_gmsh_2d_reader.cc
    "fun_model;fun_discretization"
  )
  add_gtest(cartesian_unstruct_2d_builder_test
    impl/builder/cartesian/test_cartesian_unstruct_2d_builder.cc
    "fun_model;fun_discretization"
  )
endif()

add_gtest(gmsh_3d_reader_test
  impl/builder/gmsh/test_gmsh_3d_reader.cc
  "fun_model;fun_solver;fun_discretization"
)
