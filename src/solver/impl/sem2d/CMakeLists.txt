# 2-D acoustic SEM solver (COMPILE_SEM_2D): quadrilaterals in the (x, z) plane on a
# 2-D ModelUnstruct. Header-only class, explicitly instantiated once per
# (order, model location) in generated sources.
target_compile_definitions(fun_solver PUBLIC COMPILE_SEM_2D)

target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

foreach(ORDER RANGE 1 ${TARGET_SOLVER_ACOUSTIC_ORDER})
  foreach(MODEL_TYPE onelements onnodes)
    if(MODEL_TYPE STREQUAL "onnodes")
      set(MODEL_ON_NODES "true")
    else()
      set(MODEL_ON_NODES "false")
    endif()
    set(_gen "${CMAKE_CURRENT_BINARY_DIR}/SEM2DQ${ORDER}_unstruct_makutu_solver_acoustic_model${MODEL_TYPE}.cpp")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/templates/sem_solver_2d_acoustic.cpp.in ${_gen} @ONLY)
    target_sources(fun_solver PRIVATE ${_gen})
  endforeach()
endforeach()
unset(_gen)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
