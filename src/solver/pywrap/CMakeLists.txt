if(ENABLE_PYWRAP)
  pybind11_add_module(pyproxys src/bindings.cpp)

  target_include_directories(pyproxys
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(pyproxys
    PUBLIC
      proxy_solver
  )

  configure_python_module(pyproxys)
endif()