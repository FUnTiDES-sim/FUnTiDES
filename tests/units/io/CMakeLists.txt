# test_adios2_io_controller.cc is built in every configuration: it uses only
# the public API and skips its tests when the ADIOS2 backend is not built in.
set(IO_TEST_SOURCES
  impl/test_posix_io_controller.cc
  impl/test_adios2_io_controller.cc
)

add_gtest_group(io_tests
  SOURCES ${IO_TEST_SOURCES}
  EXTRA_LINKS fun_io
)
