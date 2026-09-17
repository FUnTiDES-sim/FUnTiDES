set(PARALLEL_TEST_SOURCES
  test_topology_factory.cc
)

add_gtest_group(parallel_tests
  SOURCES ${PARALLEL_TEST_SOURCES}
  EXTRA_LINKS fun_parallel
)
