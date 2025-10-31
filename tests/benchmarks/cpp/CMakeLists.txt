add_benchmark(bench_solver_fe_cartesian_struct_acoustic
 SOURCES
   src/bench_solver_fe_cartesian_struct_acoustic.cc
 INCLUDES
   ${CMAKE_CURRENT_SOURCE_DIR}/include
 LIBS
   proxy_solver
   proxy_model_builder_cartesian
   proxy_model_struct
   discretization
   proxy_utils
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
   proxy_solver
   proxy_model_builder_cartesian
   proxy_model_struct
   discretization
   proxy_utils
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
   proxy_solver
   proxy_model_builder_cartesian
   proxy_model_unstruct
   discretization
   proxy_utils
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
   proxy_solver
   proxy_model_builder_cartesian
   proxy_model_unstruct
   discretization
   proxy_utils
 THREADS
   1 2 4 8 16 32 64
 LABELS
   unstruct
   solver
)