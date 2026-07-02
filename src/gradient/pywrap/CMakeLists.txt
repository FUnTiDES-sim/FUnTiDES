if(ENABLE_PYWRAP)
  pybind11_add_module(gradient src/bindings.cpp)

  target_include_directories(gradient
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(gradient
    PUBLIC
      fun_gradient
  )

  configure_python_module(gradient)
endif()
