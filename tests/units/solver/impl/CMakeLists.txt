link_libraries(GTest::gtest GTest::gtest_main)

add_gtest(wavefield_acoustic_test
  test_wavefield_acoustic.cc
  proxy_solver
)

add_gtest(wavefield_elastic_test
  test_wavefield_elastic.cc
  proxy_solver
)

add_gtest(rhs_acoustic_test
  test_rhs_acoustic.cc
  proxy_solver
)

add_gtest(rhs_elastic_test
  test_rhs_elastic.cc
  proxy_solver
)

if(USE_MPI)
add_executable(test_mpi_partitioning test_mpi_partitioning.cc)
target_link_libraries(test_mpi_partitioning
    PRIVATE
    proxy_model_builder_cartesian
    MPI::MPI_CXX
    GTest::gtest
    GTest::gtest_main
)

add_test(NAME test_mpi_partitioning COMMAND test_mpi_partitioning)
endif()
