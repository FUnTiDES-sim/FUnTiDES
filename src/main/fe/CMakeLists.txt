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
      fun_discretization
      fun_solver
      fun_model
      fun_core
      fun_parallel
      fun_acquisition
      adios2::adios2
      ${STDCXXFS_LIB}
  )

install(TARGETS funtides-sem
        RUNTIME DESTINATION bin)

endif()
