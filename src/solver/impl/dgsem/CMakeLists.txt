# Coupled DG-SEM solver variant (composite: reuses the SEM and DG solvers). common/
# holds the coupled class, data and implementation headers; acoustic/ the
# per-physics types and explicit instantiations.
target_compile_definitions(fun_solver PUBLIC
    MAX_DG_SEM_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

if(COMPILE_DG_SEM)
  target_compile_definitions(fun_solver PUBLIC COMPILE_DG_SEM)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

add_subdirectory(common)
add_subdirectory(acoustic)
