#-------------------------------------------------------------------
# Google benchmark configuration
#-------------------------------------------------------------------
find_package(benchmark QUIET)
if(NOT benchmark_FOUND)
  message(FATAL_ERROR
    "FUNTIDES_ENABLE_TESTS=ON requires Google Benchmark, which was not found. "
    "Install the TPL without --skip-tests (and add its prefix to CMAKE_PREFIX_PATH), "
    "or configure with -DFUNTIDES_ENABLE_TESTS=OFF.")
endif()

# Set output directory for current benchmark results (available globally)
set(BENCHMARK_RESULTS_DIR ${CMAKE_BINARY_DIR}/Benchmarking/cpp CACHE PATH "Directory for benchmark results")
file(MAKE_DIRECTORY ${BENCHMARK_RESULTS_DIR})
