#-------------------------------------------------------------------
# Google test/benchmark related functions
#-------------------------------------------------------------------

# Helper function to create a benchmark executable
# The first argument is the target name, the rest are sources, includes and libraries
# Example :
# add_benchmark(bench_solver_struct
#   SOURCES
#     src/my_bench.cc
#   INCLUDES
#     ${CMAKE_CURRENT_SOURCE_DIR}/include
#   LIBS
#     proxy_solver
#     proxy_model_builder_cartesian
#     proxy_model_struct
#     proxy_model_unstruct
#     discretization
#     proxy_utils
#   LABELS 
#     solver
#     struct
# )
function(add_benchmark name)
  set(options)      # cmake flags for verbose, etc. (optional)
  set(oneValueArgs) # name
  set(multiValueArgs SOURCES INCLUDES LIBS LABELS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  add_executable(${name} ${ARG_SOURCES})

  target_link_libraries(${name}
    PRIVATE
      ${ARG_LIBS}
      benchmark::benchmark
  )

  # Add include directories if provided
  if(ARG_INCLUDES)
    target_include_directories(${name}
      PRIVATE
        ${ARG_INCLUDES}
    )
  endif()

  target_link_kokkos_if_enabled(${name})

  set_target_properties(${name}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
  )

  add_test(
    NAME ${name}
    COMMAND $<TARGET_FILE:${name}> 
      --benchmark_out=${BENCHMARK_RESULTS_DIR}/${name}_latest.json
      --benchmark_out_format=json
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