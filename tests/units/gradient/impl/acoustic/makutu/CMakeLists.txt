# Test Struct
add_gtest(differentiator_acoustic_makutu_struct_test
  test_differentiator_acoustic_makutu_struct.cc
  funtides_gradient
)

target_compile_definitions(differentiator_acoustic_makutu_struct_test PRIVATE
    MAX_DIFFERENTIATOR_ACOUSTIC_ORDER=${TARGET_DIFFERENTIATOR_ACOUSTIC_ORDER}
)

# Test Unstruct
add_gtest(differentiator_acoustic_makutu_unstruct_test
  test_differentiator_acoustic_makutu_unstruct.cc
  funtides_gradient
)

target_compile_definitions(differentiator_acoustic_makutu_unstruct_test PRIVATE
    MAX_DIFFERENTIATOR_ACOUSTIC_ORDER=${TARGET_DIFFERENTIATOR_ACOUSTIC_ORDER}
)
