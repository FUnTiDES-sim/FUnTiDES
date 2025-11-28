if(COMPILE_FD)

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

  add_executable(fdproxy
                 src/main.cc
                 src/fd_proxy.cc
                )

  target_include_directories(fdproxy PRIVATE
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${CMAKE_SOURCE_DIR}/external/cxxopts/
                            )

  target_link_libraries(fdproxy
    PRIVATE
      discretization
      proxy_utils
      fd_grid
      fd_solver
      fd_discretization
  )

  install(TARGETS fdproxy
        RUNTIME DESTINATION bin)
endif()
