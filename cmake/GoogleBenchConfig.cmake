#-------------------------------------------------------------------
# Google benchmark configuration
#-------------------------------------------------------------------
find_package(benchmark REQUIRED)

# Set output directory for current benchmark results (available globally)
set(BENCHMARK_RESULTS_DIR ${CMAKE_BINARY_DIR}/Benchmarking/cpp CACHE PATH "Directory for benchmark results")
file(MAKE_DIRECTORY ${BENCHMARK_RESULTS_DIR})
