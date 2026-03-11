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

add_gtest(attenuation_test
  test_attenuation.cc
  proxy_solver
)
target_link_libraries(attenuation_test PRIVATE proxy_model_builder_cartesian)

add_gtest(attenuation_receiver_test
  test_attenuation_receiver.cc
  proxy_solver
)
target_link_libraries(attenuation_receiver_test PRIVATE proxy_model_builder_cartesian)