configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/bench_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/bench_config.h
  @ONLY
)

add_benchmark(bench_solver_fe_cartesian_struct_acoustic
 SOURCES
   src/bench_solver_fe_cartesian_struct_acoustic.cc
 INCLUDES
   ${CMAKE_CURRENT_SOURCE_DIR}/include
 LIBS
   fun_solver
   fun_model
   fun_discretization
   fun_utils
 THREADS
   1 2 4 8 16 32 64
 LABELS
   struct
   solver
)

add_benchmark(bench_solver_fe_cartesian_struct_elastic
 SOURCES
   src/bench_solver_fe_cartesian_struct_elastic.cc
 INCLUDES
   ${CMAKE_CURRENT_SOURCE_DIR}/include
 LIBS
   fun_solver
   fun_model
   fun_discretization
   fun_utils
 THREADS
   1 2 4 8 16 32 64
 LABELS
   struct
   solver
)


add_benchmark(bench_solver_fe_cartesian_unstruct_acoustic
 SOURCES
   src/bench_solver_fe_cartesian_unstruct_acoustic.cc
 INCLUDES
   ${CMAKE_CURRENT_SOURCE_DIR}/include
 LIBS
   fun_solver
   fun_model
   fun_discretization
   fun_utils
 THREADS
   1 2 4 8 16 32 64
 LABELS
   unstruct
   solver
)

add_benchmark(bench_solver_fe_cartesian_unstruct_elastic
 SOURCES
   src/bench_solver_fe_cartesian_unstruct_elastic.cc
 INCLUDES
   ${CMAKE_CURRENT_SOURCE_DIR}/include
 LIBS
   fun_solver
   fun_model
   fun_discretization
   fun_utils
 THREADS
   1 2 4 8 16 32 64
 LABELS
   unstruct
   solver
)

add_benchmark(bench_solver_acoustic
 SOURCES
   src/bench_acoustic_kernel.cc
 INCLUDES
   ${CMAKE_CURRENT_SOURCE_DIR}/include
 LIBS
   fun_solver
   fun_model
   fun_discretization
   fun_utils
 THREADS
   1
 LABELS
   unstruct
   solver
)

# Standalone C-PML acoustic SEM proxy harness (not a Google Benchmark binary:
# the proxy runner invokes it as <binary> <param_file>). Built unconditionally
# like the other struct benchmarks.
add_executable(bench_solver_fe_cartesian_struct_acoustic_pml
  src/bench_solver_fe_cartesian_struct_acoustic_pml.cc
)
target_link_libraries(bench_solver_fe_cartesian_struct_acoustic_pml
  PRIVATE
    fun_solver
    fun_model
    fun_discretization
    fun_utils
)
target_include_directories(bench_solver_fe_cartesian_struct_acoustic_pml
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_kokkos_if_enabled(bench_solver_fe_cartesian_struct_acoustic_pml)
set_target_properties(bench_solver_fe_cartesian_struct_acoustic_pml
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
)
