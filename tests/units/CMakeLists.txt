# Set google test main file with kokkos initialization/finalization
set(TEST_MAIN ${CMAKE_SOURCE_DIR}/tests/units/main.cpp)

add_subdirectory(model)
add_subdirectory(utils)