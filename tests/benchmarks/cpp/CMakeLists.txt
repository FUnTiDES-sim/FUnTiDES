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
   fun_acquisition
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
   fun_acquisition
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
   fun_acquisition
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
   fun_acquisition
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
   fun_acquisition
 THREADS
   1
 LABELS
   unstruct
   solver
)

# ---------------------------------------------------------------------------
# MIMIR proxy harnesses: standalone executables (no google-benchmark), driven
# by the proxy server. Each prints a metrics block (time_s=, finite=, ...) and
# writes the final field as raw float64 to the file passed as argv[1]. The
# build stays non-O3: the SEM solver has a latent use-after-free that -O3
# turns into a crash (see scripts/build_sem_proxy.sh).
# ---------------------------------------------------------------------------
add_executable(proxy_sem_acoustic src/proxy_sem_acoustic.cc)
target_include_directories(proxy_sem_acoustic
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
)
target_link_libraries(proxy_sem_acoustic
  PRIVATE
    fun_solver
    fun_model
    fun_discretization
    fun_acquisition
)
target_link_kokkos_if_enabled(proxy_sem_acoustic)
set_target_properties(proxy_sem_acoustic
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
)

add_executable(proxy_sem_elastic src/proxy_sem_elastic.cc)
target_include_directories(proxy_sem_elastic
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
)
target_link_libraries(proxy_sem_elastic
  PRIVATE
    fun_solver
    fun_model
    fun_discretization
    fun_acquisition
)
target_link_kokkos_if_enabled(proxy_sem_elastic)
set_target_properties(proxy_sem_elastic
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
)
