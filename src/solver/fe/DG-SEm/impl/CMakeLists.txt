target_compile_definitions(fun_solver PUBLIC
    MAX_DG_SEM_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
)

if(COMPILE_DG_SEM)
  target_compile_definitions(fun_solver PUBLIC COMPILE_DG_SEM)
endif()

add_subdirectory(acoustic)
add_subdirectory(common)
