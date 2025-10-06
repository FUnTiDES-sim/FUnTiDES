#-------------------------------------------------------------------
# Google test/benchmark configuration
#-------------------------------------------------------------------

# Googletest is bundled in benchmark
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/googletest/CMakeLists.txt")
    message(FATAL_ERROR "googletest submodule not found. Run 'git submodule update --init --recursive'")
else()
    add_subdirectory(${CMAKE_SOURCE_DIR}/external/googletest)
    include(GoogleTest)
    enable_testing()
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/external/benchmark/CMakeLists.txt")
    set(GOOGLETEST_PATH "${CMAKE_SOURCE_DIR}/external/googletest")
    set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
    set(BENCHMARK_ENABLE_TESTING OFF)
    add_subdirectory(${CMAKE_SOURCE_DIR}/external/benchmark)
    # Set output directory for current benchmark results (available globally)
    set(BENCHMARK_RESULTS_DIR ${CMAKE_BINARY_DIR}/Benchmarking CACHE PATH "Directory for benchmark results")
    file(MAKE_DIRECTORY ${BENCHMARK_RESULTS_DIR})
else()
    message(FATAL_ERROR "benchmark submodule not found. Run 'git submodule update --init --recursive'")
endif()