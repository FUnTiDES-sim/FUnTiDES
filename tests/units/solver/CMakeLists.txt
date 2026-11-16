link_libraries(GTest::gtest GTest::gtest_main)

# The tree mirrors src/solver: the shared per-physics types under
# impl/common/<physics>, one subtree per method under impl/. Method-specific
# sources compile only under their COMPILE_* option here, and carry the
# matching #ifdef guard in the source.
set(SOLVER_TEST_SOURCES
  impl/common/acoustic/test_rhs_acoustic.cc
  impl/common/elastic/test_rhs_elastic.cc
  impl/common/acoustoelastic/test_rhs_acoustoelastic.cc
)

if(COMPILE_SEM)
  list(APPEND SOLVER_TEST_SOURCES
    impl/common/test_solver_factory.cc
    impl/sem/common/test_physics_traits.cc
    impl/sem/common/test_solver_data_acoustic.cc
    impl/sem/common/test_solver_data_elastic.cc
    impl/sem/common/test_attenuation.cc
    impl/sem/common/test_boundary_synchronizer.cc
    impl/sem/common/test_sem_solver_zdeformed.cc
    impl/sem/acoustic/test_wavefield_acoustic.cc
    impl/sem/acoustic/test_sem_solver_acoustic.cc
    impl/sem/elastic/test_wavefield_elastic.cc
    impl/sem/elastic/test_sem_solver_elastic.cc
    impl/sem/elastic/test_elastic_flux_symmetry.cc
    impl/sem/elastic/test_elastic_flux_tti_compact.cc
    impl/sem/acoustoelastic/test_wavefield_acoustoelastic.cc
    impl/sem/acoustoelastic/test_sem_solver_acoustoelastic.cc
  )
endif()

if(COMPILE_DG)
  list(APPEND SOLVER_TEST_SOURCES
    impl/dg/test_dg_solver_factory.cc
    impl/dg/test_dg_penalty.cc
    impl/dg/test_dg_solver_data.cc
  )
endif()

if(COMPILE_DG_SEM)
  list(APPEND SOLVER_TEST_SOURCES
    impl/dgsem/test_dg_sem_solver_factory.cc
    impl/dgsem/test_dgsem_solver_data.cc
    impl/dgsem/test_dgsem_solver_acoustic.cc
  )
endif()

if(COMPILE_DG_PADAPTIVE)
  list(APPEND SOLVER_TEST_SOURCES
    impl/dg_padaptive/test_dg_padaptive_solver_factory.cc
    impl/dg_padaptive/test_dg_padaptive_solver_data.cc
    impl/dg_padaptive/test_dg_padaptive_solver_acoustic.cc
  )
endif()

add_gtest_group(solver_tests
  SOURCES ${SOLVER_TEST_SOURCES}
  EXTRA_LINKS fun_solver fun_model
)

# The z-deformed team kernels need one thread per element node and fall back to
# simpler paths below that, so their tests run once more with enough threads.
# OMP_THREAD_LIMIT is raised too: below the Kokkos thread count, team barriers
# wait for threads that are never created.
add_test(NAME solver_tests_zdeformed_team
  COMMAND solver_tests --gtest_filter=SemSolverZDeformed.*
)
set_tests_properties(solver_tests_zdeformed_team PROPERTIES
  ENVIRONMENT "OMP_NUM_THREADS=12;OMP_THREAD_LIMIT=12;KOKKOS_NUM_THREADS=12;OMP_PROC_BIND=false;OMP_WAIT_POLICY=passive"
  TIMEOUT 120
)

# ============================================================================
# One DGsolver acoustic test binary per polynomial order: each order is a
# distinct explicit template instantiation of DGsolver (see
# dg_solver_impl.h), so splitting per order (like the Qk_Hexahedron
# discretization tests) is what actually exercises order > 1.
# ============================================================================
if(COMPILE_DG)
  foreach(ORDER RANGE 1 ${TARGET_SOLVER_ACOUSTIC_ORDER})
    set(_gen_file "${CMAKE_CURRENT_BINARY_DIR}/test_dg_solver_acoustic_Q${ORDER}.cpp")

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/impl/dg/templates/test_dg_solver_acoustic_order.cpp.in"
      "${_gen_file}"
      @ONLY
    )

    set(_target "dg_solver_acoustic_Q${ORDER}_test")
    add_gtest(${_target}
      ${_gen_file}
      "fun_solver;fun_model"
    )

    target_include_directories(${_target} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/impl/dg
    )
  endforeach()

  unset(_gen_file)
  unset(_target)
endif()

if(USE_MPI)
add_executable(test_mpi_partitioning test_mpi_partitioning.cc)
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
# impl/sem/common/test_attenuation_receiver.cc,
# impl/sem/common/test_attenuation_reference.cc
