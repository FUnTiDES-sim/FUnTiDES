if(COMPILE_FD)
  add_executable(fdproxy
                 src/main.cc
                 src/fdtd_proxy.cc
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
      fdtd_solver
      fdtd_discretization
  )

  target_link_kokkos_if_enabled(fdproxy)
endif()
