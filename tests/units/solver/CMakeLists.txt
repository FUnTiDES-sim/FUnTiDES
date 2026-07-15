link_libraries(GTest::gtest GTest::gtest_main)

set(SOLVER_TEST_SOURCES
  impl/test_wavefield_acoustic.cc
  impl/test_wavefield_elastic.cc
  impl/test_rhs_acoustic.cc
  impl/test_rhs_elastic.cc
  impl/test_wavefield_acoustoelastic.cc
  impl/test_rhs_acoustoelastic.cc
  impl/test_sem_solver_acoustoelastic.cc
  impl/test_physics_traits.cc
  impl/test_solver_data_acoustic.cc
  impl/test_solver_data_elastic.cc
  impl/test_solver_factory.cc
  impl/test_attenuation.cc
  impl/test_sem_solver_acoustic.cc
  impl/test_sem_solver_elastic.cc
  impl/test_boundary_synchronizer.cc
)

if(COMPILE_DG)
  list(APPEND SOLVER_TEST_SOURCES
    impl/test_dg_solver_factory.cc
    dg/test_dg_penalty.cc
    dg/test_dg_solver_data.cc
    dg/test_dg_solver_acoustic.cc
  )
endif()

add_gtest_group(solver_tests
  SOURCES ${SOLVER_TEST_SOURCES}
  EXTRA_LINKS fun_solver fun_model
)

if(USE_MPI)
add_executable(test_mpi_partitioning impl/test_mpi_partitioning.cc)
target_link_libraries(test_mpi_partitioning
    PRIVATE
    fun_model
    MPI::MPI_CXX
    GTest::gtest
    GTest::gtest_main
)
add_test(NAME test_mpi_partitioning COMMAND test_mpi_partitioning)
endif()

# Disabled — kept for reference, not wired into any target:
# impl/test_attenuation_receiver.cc, impl/test_attenuation_reference.cc
