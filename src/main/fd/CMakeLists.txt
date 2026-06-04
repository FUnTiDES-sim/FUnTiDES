if(COMPILE_FD)

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

  add_executable(funtides-fd
                 src/main.cc
                 src/fd_proxy.cc
                )

  target_include_directories(funtides-fd PRIVATE
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${CMAKE_SOURCE_DIR}/external/cxxopts/
                            )

  target_link_libraries(funtides-fd
    PRIVATE
      fun_discretization
      fun_utils
      fd_grid
      fd_solver
      fd_discretization
  )

  install(TARGETS funtides-fd
        RUNTIME DESTINATION bin)
endif()
