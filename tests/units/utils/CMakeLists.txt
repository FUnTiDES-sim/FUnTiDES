add_gtest(common_macros_find_max_1d_test
  test_common_macros_find_max_1d.cc
  fun_utils
)

add_gtest(das_receiver_test
  test_das_receiver.cc
  discretization
)

add_gtest(evaluate_ricker_test
  test_evaluate_ricker.cc
  fun_utils
)

add_gtest(source_and_receiver_utils_test
  test_source_and_receiver_utils.cc
  discretization
)

add_gtest(args_parse_test
  test_args_parse.cc
  fun_utils
)