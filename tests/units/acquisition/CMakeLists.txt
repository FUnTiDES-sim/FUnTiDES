set(ACQUISITION_TEST_SOURCES
  test_das_receiver.cc
  test_evaluate_ricker.cc
  test_source_and_receiver_utils.cc
)

add_gtest_group(acquisition_tests
  SOURCES ${ACQUISITION_TEST_SOURCES}
  EXTRA_LINKS fun_acquisition
)
