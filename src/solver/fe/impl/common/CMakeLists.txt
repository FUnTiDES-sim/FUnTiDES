target_sources(fun_solver PRIVATE 
     src/solver_factory.cc
)

target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
