link_libraries(GTest::gtest GTest::gtest_main)

add_gtest(dg-sem_solver_data_test
  test_dg-sem_solver_data.cc
  fun_solver
)

add_gtest(dg-sem_solver_acoustic_test
  test_dg-sem_solver_acoustic.cc
  fun_solver
)
target_link_libraries(dg-sem_solver_acoustic_test PRIVATE fun_model)
