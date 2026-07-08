add_gtest(mesh_struct_test
  test_struct_mesh.cc
  fun_model
)

# Restrict the polynomial orders exercised by the typed tests to those actually
# compiled on this machine (higher orders are not built here).
target_compile_definitions(mesh_struct_test PRIVATE
    MAX_SOLVER_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

add_gtest(gllpoints_test
  test_gllpoints.cc
  fun_model
)
