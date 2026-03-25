add_gtest(gradient_acoustic_test
  test_gradient_acoustic.cc
  funtides_gradient
)

add_gtest(wavefield_view_forward_acoustic_test
  test_wavefield_view_forward_acoustic.cc
  funtides_gradient
)

add_gtest(wavefield_view_backward_acoustic_test
  test_wavefield_view_backward_acoustic.cc
  funtides_gradient
)

add_subdirectory(makutu)
