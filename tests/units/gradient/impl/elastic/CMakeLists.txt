add_gtest(gradient_data_elastic_test
  test_gradient_data_elastic.cc
  fun_gradient
)

add_gtest(gradient_elastic_test
  test_gradient_elastic.cc
  fun_gradient
)

add_gtest(wavefield_view_forward_elastic_test
  test_wavefield_view_forward_elastic.cc
  fun_gradient
)

add_gtest(wavefield_view_backward_elastic_test
  test_wavefield_view_backward_elastic.cc
  fun_gradient
)

add_subdirectory(makutu)
