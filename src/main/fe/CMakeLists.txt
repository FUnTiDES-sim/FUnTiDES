if(COMPILE_SEM)

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

  add_executable(funtides-sem
    src/main.cc
    src/sem_proxy.cc
  )

  target_include_directories(funtides-sem
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${CMAKE_SOURCE_DIR}/external/cxxopts
  )

  target_link_libraries(funtides-sem
    PRIVATE
      discretization
      proxy_solver
      proxy_model_builder_cartesian
      proxy_model_struct
      proxy_model_unstruct
      proxy_utils
      adios2::adios2
      ${STDCXXFS_LIB}
  )

install(TARGETS funtides-sem
        RUNTIME DESTINATION bin)

endif()
