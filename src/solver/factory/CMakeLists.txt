# Construction entry point shared by every solver method. solver_factory.cc
# was lifted out of the SEM impl/ tree (review 2026-10-19): the factory is
# method-agnostic, so it lives at the common level beside api/, and its
# header joins the target's public include path.
target_sources(fun_solver
  PRIVATE
    src/solver_factory.cc
)

target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
