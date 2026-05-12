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

add_gtest(solver_data_test
  test_solver_data.cc
  proxy_solver
)

add_gtest(solver_factory_test
  test_solver_factory.cc
  proxy_solver
)

add_gtest(attenuation_test
  test_attenuation.cc
  proxy_solver
)
target_link_libraries(attenuation_test PRIVATE proxy_model_builder_cartesian)

add_gtest(sem_solver_acoustic_test
  test_sem_solver_acoustic.cc
  proxy_solver
)
target_link_libraries(sem_solver_acoustic_test PRIVATE proxy_model_builder_cartesian)

add_gtest(sem_solver_elastic_test
  test_sem_solver_elastic.cc
  proxy_solver
)
target_link_libraries(sem_solver_elastic_test PRIVATE proxy_model_builder_cartesian)

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


# target_link_libraries(attenuation_test PRIVATE proxy_model_builder_cartesian)

# add_gtest(attenuation_receiver_test
#   test_attenuation_receiver.cc
#   proxy_solver
# )
# target_link_libraries(attenuation_receiver_test PRIVATE proxy_model_builder_cartesian)

# add_gtest(attenuation_reference_test
#   test_attenuation_reference.cc
#   proxy_solver
# )
# target_link_libraries(attenuation_reference_test PRIVATE proxy_model_builder_cartesian)
