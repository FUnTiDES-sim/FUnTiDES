if (ENABLE_PYWRAP)

  if(COMPILE_SEM)

    file(COPY ${CMAKE_SOURCE_DIR}/scripts/pyscript/ DESTINATION ${CMAKE_BINARY_DIR}/proxy)
    pybind11_add_module(proxy_pysem src/bindings.cpp)

    target_include_directories(proxy_pysem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

    target_link_libraries(proxy_pysem
      PRIVATE
        proxy_solver
        proxy_model_builder_cartesian
        proxy_model_struct
        proxy_model_unstruct
        proxy_utils
    )

    configure_python_module(proxy_pysem)

  endif()

endif()