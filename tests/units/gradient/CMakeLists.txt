set(GRADIENT_TEST_SOURCES
  impl/common/test_differentiator_factory.cc
  impl/acoustic/test_gradient_acoustic.cc
  impl/acoustic/test_gradient_data_acoustic.cc
  impl/acoustic/test_wavefield_view_forward_acoustic.cc
  impl/acoustic/test_wavefield_view_backward_acoustic.cc
  impl/acoustic/makutu/test_differentiator_acoustic_makutu_struct.cc
  impl/acoustic/makutu/test_differentiator_acoustic_makutu_unstruct.cc
  impl/elastic/test_gradient_data_elastic.cc
  impl/elastic/test_gradient_elastic.cc
  impl/elastic/test_wavefield_view_forward_elastic.cc
  impl/elastic/test_wavefield_view_backward_elastic.cc
  impl/elastic/makutu/test_differentiator_elastic_makutu_struct.cc
  impl/elastic/makutu/test_differentiator_elastic_makutu_unstruct.cc
)

add_gtest_group(gradient_tests
  SOURCES ${GRADIENT_TEST_SOURCES}
  EXTRA_LINKS fun_gradient
)

target_compile_definitions(gradient_tests PRIVATE
    MAX_DIFFERENTIATOR_ACOUSTIC_ORDER=${TARGET_DIFFERENTIATOR_ACOUSTIC_ORDER}
    MAX_DIFFERENTIATOR_ELASTIC_ORDER=${TARGET_DIFFERENTIATOR_ELASTIC_ORDER}
)
