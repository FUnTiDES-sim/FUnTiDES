if(COMPILE_SEM)
  add_executable(semproxy
    src/main.cpp
    src/SEMproxy.cpp
  )

  target_include_directories(semproxy
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${CMAKE_SOURCE_DIR}/external/cxxopts/
  )

  target_link_libraries(semproxy
    PRIVATE
      discretization
      solver
      proxy_model_builder_cartesian
      proxy_model_struct
      proxy_model_unstruct
      proxy_utils
  )

  if (ENABLE_PYWRAP)
    file(COPY ${CMAKE_SOURCE_DIR}/scripts/pyscript/ DESTINATION ${CMAKE_BINARY_DIR}/proxy)
    pybind11_add_module(pysem src/pysem.cpp src/SEMproxy.cpp)

    target_include_directories(pysem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

    target_link_libraries(pysem
      PRIVATE
        solver
        proxy_model_builder_cartesian
        proxy_model_struct
        proxy_model_unstruct
        proxy_utils
    )

    configure_python_module(pysem)
  endif()
endif()


################################
# Setup for FD_EXECUTABLE
################################
if(COMPILE_FD)
  add_executable(fdproxy src/FDTDmain.cpp)
  # target_include_directories(fd_executable PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/utilsFDTD)
  target_include_directories(fdproxy PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../discretization
  )

  target_link_libraries(fdproxy
    PRIVATE
      discretization
      proxy_utils
  )

  target_link_kokkos_if_enabled(fdproxy)
endif()
