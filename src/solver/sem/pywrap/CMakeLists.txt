if(ENABLE_PYWRAP)
  pybind11_add_module(solver src/bindings.cpp)

  target_include_directories(solver
    PRIVATE
      include
      ${CMAKE_CURRENT_SOURCE_DIR}/../dg/pywrap/include
      ${CMAKE_CURRENT_SOURCE_DIR}/../dgsem/pywrap/include
      ${CMAKE_CURRENT_SOURCE_DIR}/../dg_padaptive/pywrap/include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(solver
    PUBLIC
      fun_solver
  )

  configure_python_module(solver)
endif()
