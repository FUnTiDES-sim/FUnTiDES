if(ENABLE_PYWRAP)
  pybind11_add_module(proxy_pysolver src/bindings.cpp)

  target_include_directories(proxy_pysolver
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(proxy_pysolver
    PUBLIC
      proxy_solver
  )

  configure_python_module(proxy_pysolver)
endif()