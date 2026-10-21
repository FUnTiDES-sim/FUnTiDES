# DG solver variant. common/ holds the class, data and penalty headers with
# the implementation header; acoustic/ holds the per-physics types and the
# explicit instantiations. All variants compile into fun_solver (declared in
# impl/CMakeLists.txt), so this directory only contributes definitions and sources.
target_compile_definitions(fun_solver PUBLIC
    MAX_DG_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

if(COMPILE_DG)
  target_compile_definitions(fun_solver PUBLIC COMPILE_DG)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

add_subdirectory(common)
add_subdirectory(acoustic)
