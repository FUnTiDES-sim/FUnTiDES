if(ENABLE_PYWRAP)
  pybind11_add_module(io src/bindings.cpp)

  target_include_directories(io
    PRIVATE
      include
  )

  # fun_io carries the ADIOS2 backend when FUNTIDES_ENABLE_ADIOS2=ON; the
  # module itself never includes ADIOS2 headers.
  target_link_libraries(io PUBLIC fun_io)

  configure_python_module(io)
endif()
