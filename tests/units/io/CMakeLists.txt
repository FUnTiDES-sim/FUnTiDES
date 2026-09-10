set(IO_TEST_SOURCES
  impl/test_posix_io_controller.cc
)

add_gtest_group(io_tests
  SOURCES ${IO_TEST_SOURCES}
  EXTRA_LINKS funtides_io
)
