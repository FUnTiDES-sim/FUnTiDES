target_sources(proxy_solver PRIVATE 
     src/solver_factory.cc
)

target_include_directories(proxy_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
