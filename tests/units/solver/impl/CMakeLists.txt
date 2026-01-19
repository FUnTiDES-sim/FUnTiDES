add_gtest(wavefield_acoustic_test
  test_wavefield_acoustic.cc
  proxy_solver
)

add_gtest(wavefield_elastic_test
  test_wavefield_elastic.cc
  proxy_solver
)

add_gtest(rhs_acoustic_test
  test_rhs_acoustic.cc
  proxy_solver
)

add_gtest(rhs_elastic_test
  test_rhs_elastic.cc
  proxy_solver
)