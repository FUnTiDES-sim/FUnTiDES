set(CORE_TEST_SOURCES
  test_common_macros_find_max_1d.cc
)

add_gtest_group(core_tests
  SOURCES ${CORE_TEST_SOURCES}
  EXTRA_LINKS fun_core
)
