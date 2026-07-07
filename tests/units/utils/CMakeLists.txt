set(UTILS_TEST_SOURCES
  test_common_macros_find_max_1d.cc
  test_das_receiver.cc
  test_evaluate_ricker.cc
  test_source_and_receiver_utils.cc
  test_args_parse.cc
)

add_gtest_group(utils_tests
  SOURCES ${UTILS_TEST_SOURCES}
  EXTRA_LINKS fun_utils fun_discretization
)
