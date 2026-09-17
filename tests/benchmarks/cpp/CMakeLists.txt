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

# Standalone C-PML elastic SEM proxy harness (not a Google Benchmark binary:
# the proxy runner invokes it as <binary> <param_file>). Built unconditionally
# like the other struct benchmarks.
#
# The harness instantiates SEMsolverElastic<1, MAKUTU, ...> in its own TU, so
# it does NOT need the pre-instantiated fun_solver library (104 generated TU
# that rebuild on every solver-header touch). It links only fun_model (real
# symbols) plus the header-only INTERFACE targets for the include paths, and
# adds the solver-impl include dirs explicitly. This keeps each ratchet build
# to a single TU instead of the whole fun_solver library.
add_executable(bench_solver_fe_cartesian_struct_elastic_pml
  src/bench_solver_fe_cartesian_struct_elastic_pml.cc
)
target_link_libraries(bench_solver_fe_cartesian_struct_elastic_pml
  PRIVATE
    fun_model
    fun_discretization_impl
    fun_solver_api
)
target_include_directories(bench_solver_fe_cartesian_struct_elastic_pml
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/solver/fe/impl/common/include
    ${CMAKE_SOURCE_DIR}/src/solver/fe/impl/acoustic/include
    ${CMAKE_SOURCE_DIR}/src/solver/fe/impl/elastic/include
)
target_link_kokkos_if_enabled(bench_solver_fe_cartesian_struct_elastic_pml)
set_target_properties(bench_solver_fe_cartesian_struct_elastic_pml
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
)
