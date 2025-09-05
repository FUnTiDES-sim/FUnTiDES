if(ENABLE_PYWRAP)
  pybind11_add_module(solver src/bindings.cpp)

  target_include_directories(solver
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(solver
    PUBLIC
      proxy_solver
  )

  configure_python_module(solver)
endif()