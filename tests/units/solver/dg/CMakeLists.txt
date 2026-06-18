link_libraries(GTest::gtest GTest::gtest_main)

add_gtest(dg_penalty_test
  test_dg_penalty.cc
  fun_solver
)

add_gtest(dg_solver_data_test
  test_dg_solver_data.cc
  fun_solver
)

add_gtest(dg_solver_acoustic_test
  test_dg_solver_acoustic.cc
  fun_solver
)
target_link_libraries(dg_solver_acoustic_test PRIVATE fun_model)

add_gtest(dg_boundary_comm_test
  ${CMAKE_CURRENT_SOURCE_DIR}/test_dg_boundary_comm.cc
  fun_solver
)
