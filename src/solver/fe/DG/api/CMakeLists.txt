# Public class/data headers of the DG solver variant. Mirrors fe/api: the
# class header is consumed by solver_factory.cc, the data header by
# src/main/fe (sem_proxy); includes resolve by bare name via this PUBLIC
# include directory, so no consumer path changes.
target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
