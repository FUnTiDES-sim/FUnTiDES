add_gtest(gradient_acoustic_test
  test_gradient_acoustic.cc
  fun_gradient
)

add_gtest(gradient_data_acoustic_test
  test_gradient_data_acoustic.cc
  fun_gradient
)

add_gtest(wavefield_view_forward_acoustic_test
  test_wavefield_view_forward_acoustic.cc
  fun_gradient
)

add_gtest(wavefield_view_backward_acoustic_test
  test_wavefield_view_backward_acoustic.cc
  fun_gradient
)

add_subdirectory(makutu)
