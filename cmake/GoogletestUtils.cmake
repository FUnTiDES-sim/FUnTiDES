#-------------------------------------------------------------------
# Google test/benchmark related functions
#-------------------------------------------------------------------

# Helper function to create a benchmark executable and register tests
#
# Creates a benchmark executable, links necessary libraries, and registers
# CTest entries for different thread counts.
#
# Usage:
#   add_benchmark(<name>
#     SOURCES <source_files>...
#     [INCLUDES <include_dirs>...]
#     [LIBS <libraries>...]
#     [THREADS <thread_counts>...]
#     [LABELS <test_labels>...]
#   )
#
# Parameters:
#   name                - Name of the benchmark executable
#   SOURCES            - List of source files for the benchmark (required)
#   INCLUDES           - List of include directories (optional)
#   LIBS               - List of libraries to link against (optional)
#   THREADS            - List of thread counts to test with (optional, default: 1)
#                        Creates separate test entries for each count
#   LABELS             - List of custom labels for CTest (optional)
#                        'benchmark' label is always added automatically
#
# Example:
#   add_benchmark(bench_solver_struct
#     SOURCES
#       src/my_bench.cc
#     INCLUDES
#       ${CMAKE_CURRENT_SOURCE_DIR}/include
#     LIBS
#       proxy_solver
#       proxy_model_builder_cartesian
#       proxy_model_struct
#       proxy_model_unstruct
#       discretization
#       proxy_utils
#     THREADS
#       2 4 8
#     LABELS 
#       solver
#       struct
#   )
function(add_benchmark name)
  set(options)      # cmake flags for verbose, etc. (optional)
  set(oneValueArgs) # name
  set(multiValueArgs SOURCES INCLUDES LIBS THREADS LABELS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT ARG_THREADS)
    set(ARG_THREADS 1)
  endif()

  if(ARG_LABELS)
    set(base_labels "${ARG_LABELS};benchmark")
  else()
    set(base_labels "benchmark")
  endif()

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

  # Register tests for different thread counts
  foreach(t IN LISTS ARG_THREADS)
    set(test_name ${name}_t${t})
  
    add_test(
      NAME ${test_name}
      COMMAND $<TARGET_FILE:${name}>
        --kokkos-threads=${t}
        --benchmark_out=${BENCHMARK_RESULTS_DIR}/${test_name}_latest.json
        --benchmark_out_format=json
    )
    set_tests_properties(${test_name}
      PROPERTIES
        LABELS "${base_labels};t${t}"
        ENVIRONMENT "OMP_NUM_THREADS=${t};KOKKOS_NUM_THREADS=${t}"
    )
  endforeach()
endfunction()