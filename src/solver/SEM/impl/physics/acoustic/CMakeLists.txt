# Acoustic physics fragment (physics-policy migration, step 1):
# sem_solver_acoustic_kernels.h, extracted verbatim from
# impl/common/include/sem_solver_kernels.h.
target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
