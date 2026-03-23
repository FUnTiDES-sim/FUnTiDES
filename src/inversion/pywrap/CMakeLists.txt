if(ENABLE_PYWRAP)
  pybind11_add_module(inversion src/bindings.cpp)

  target_include_directories(inversion
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(inversion
    PUBLIC
      funtides_inversion
  )

  configure_python_module(inversion)
endif()
