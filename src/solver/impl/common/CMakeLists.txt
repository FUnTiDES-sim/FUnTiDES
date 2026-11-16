# Shared across methods: the construction entry point (solver_factory, the
# only place that knows every method) and, per physics, the types that every
# method consumes (RHS and physics descriptors).
target_sources(fun_solver
  PRIVATE
    src/solver_factory.cc
)

target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(acoustic)
add_subdirectory(elastic)
add_subdirectory(acoustoelastic)
