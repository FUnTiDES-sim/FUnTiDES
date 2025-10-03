#-------------------------------------------------------------------
# Google test/benchmark configuration
#-------------------------------------------------------------------

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/googletest/CMakeLists.txt")
    message(FATAL_ERROR "googletest submodule not found. Run 'git submodule update --init --recursive'")
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/external/benchmark/CMakeLists.txt")
    set(GOOGLETEST_PATH "${CMAKE_SOURCE_DIR}/external/googletest")
    set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
    set(BENCHMARK_ENABLE_TESTING ON)
    set(BENCHMARK_USE_BUNDLED_GTEST ON)
    add_subdirectory(${CMAKE_SOURCE_DIR}/external/benchmark)
else()
    message(FATAL_ERROR "benchmark submodule not found. Run 'git submodule update --init --recursive'")
endif()