target_compile_definitions(fun_solver PUBLIC
    MAX_DG_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

if(COMPILE_DG)
  target_compile_definitions(fun_solver PUBLIC COMPILE_DG)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

add_subdirectory(acoustic)
add_subdirectory(common)
