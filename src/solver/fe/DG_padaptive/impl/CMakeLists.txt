target_compile_definitions(fun_solver PUBLIC
    MAX_DG_PADAPTIVE_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

if(COMPILE_DG_PADAPTIVE)
  target_compile_definitions(fun_solver PUBLIC COMPILE_DG_PADAPTIVE)
endif()

add_subdirectory(acoustic)
add_subdirectory(common)
