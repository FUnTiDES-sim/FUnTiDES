if(ENABLE_PYWRAP)
  pybind11_add_module(model src/bindings.cpp)

  target_include_directories(model
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(model
    PUBLIC
      proxy_model_unstruct
  )

  configure_python_module(model)
endif()