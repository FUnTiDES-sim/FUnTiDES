link_libraries(GTest::gtest GTest::gtest_main)

add_gtest(dg_padaptive_solver_data_test
  test_dg_padaptive_solver_data.cc
  fun_solver
)

add_gtest(dg_padaptive_solver_acoustic_test
  test_dg_padaptive_solver_acoustic.cc
  fun_solver
)
target_link_libraries(dg_padaptive_solver_acoustic_test PRIVATE fun_model)
