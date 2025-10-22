if(ENABLE_PYWRAP)
  pybind11_add_module(model src/bindings.cpp)

  target_include_directories(model
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(model
    PUBLIC
      proxy_model_struct
      proxy_model_unstruct
      proxy_model_builder_cartesian
  )

  configure_python_module(model)
endif()