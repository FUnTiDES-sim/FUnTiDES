# Physics-generic skeleton of the SEM solver: the SEMsolver class
# declaration, the umbrella and thematic fragments (step/kernels/update/
# init/element_subset) and the shared helpers.
target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
