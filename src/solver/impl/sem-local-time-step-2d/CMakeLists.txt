# Experimental source-free two-level LF-LTS, not registered in the factory yet.
target_sources(fun_solver PRIVATE src/sem_2d_local_time_step_plan.cc)
target_include_directories(fun_solver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)

foreach(ORDER RANGE 1 ${TARGET_SOLVER_ACOUSTIC_ORDER})
  foreach(MODEL_TYPE onelements onnodes)
    if(MODEL_TYPE STREQUAL "onnodes")
      set(MODEL_ON_NODES true)
    else()
      set(MODEL_ON_NODES false)
    endif()
    set(_lts_gen "${CMAKE_CURRENT_BINARY_DIR}/SEM2DLTSQ${ORDER}_${MODEL_TYPE}.cpp")
    configure_file(templates/sem_solver_2d_local_time_step.cpp.in ${_lts_gen} @ONLY)
    target_sources(fun_solver PRIVATE ${_lts_gen})
  endforeach()
endforeach()
unset(_lts_gen)

# Validation and timing programs: tests/validation/sem-local-time-step-2d.
