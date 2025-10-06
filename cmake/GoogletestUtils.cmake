#-------------------------------------------------------------------
# Google test/benchmark related functions
#-------------------------------------------------------------------

# Helper function to create a benchmark executable
# The first argument is the target name, the rest are sources and libraries
# Example :
# add_benchmark(bench_solver_struct
#   SOURCES solver_struct_bench.cpp
#   LIBS
#     proxy_solver
#     proxy_model_builder_cartesian
#     proxy_model_struct
#     proxy_model_unstruct
#     discretization
#     proxy_utils
#   TIME 0.20
#   LABELS 
#     solver
#     struct
# )
function(add_benchmark name)
  set(options)
  set(oneValueArgs TIME)
  set(multiValueArgs SOURCES LIBS LABELS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  add_executable(${name} ${ARG_SOURCES})

  target_link_libraries(${name}
    PRIVATE
      ${ARG_LIBS}
      benchmark::benchmark
  )

  target_link_kokkos_if_enabled(${name})

  set_target_properties(${name}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
  )

  # Use provided TIME or default to 0.01
  if(ARG_TIME)
    set(benchmark_time ${ARG_TIME})
  else()
    set(benchmark_time 0.01)
  endif()

  add_test(
    NAME ${name}
    COMMAND ${name} --benchmark_min_time=${benchmark_time}
  )

  # Add labels to the benchmark test if provided, always include "benchmark"
  if(ARG_LABELS)
    set(test_labels "${ARG_LABELS};benchmark")
  else()
    set(test_labels "benchmark")
  endif()
  set_tests_properties(${name}
    PROPERTIES
      LABELS "${test_labels}"
  )
endfunction()