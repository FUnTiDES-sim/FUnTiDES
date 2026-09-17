# Set google test main file with kokkos initialization/finalization
set(TEST_MAIN ${CMAKE_SOURCE_DIR}/tests/units/main.cpp)

# Convenience target: build all tests without rebuilding libraries
add_custom_target(build_tests)

add_subdirectory(core)
add_subdirectory(model)
add_subdirectory(parallel)
add_subdirectory(acquisition)
add_subdirectory(solver)
add_subdirectory(gradient)
add_subdirectory(discretization)
