# p-adaptive DG solver variant (composite: reuses dg/ at two polynomial
# orders). common/ holds the class, data and implementation headers;
# acoustic/ the per-physics types and the paired-order instantiations.
target_compile_definitions(fun_solver PUBLIC
    MAX_DG_PADAPTIVE_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

if(COMPILE_DG_PADAPTIVE)
  target_compile_definitions(fun_solver PUBLIC COMPILE_DG_PADAPTIVE)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

add_subdirectory(common)
add_subdirectory(acoustic)
