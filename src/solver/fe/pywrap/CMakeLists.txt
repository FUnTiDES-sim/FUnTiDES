if(ENABLE_PYWRAP AND NOT COMPILE_DG)
  pybind11_add_module(solver src/bindings.cpp)

  target_include_directories(solver
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(solver
    PUBLIC
      fun_solver
  )

  configure_python_module(solver)
endif()