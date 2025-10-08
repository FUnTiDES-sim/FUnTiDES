if(COMPILE_FD)

  add_executable(fdproxy 
                 src/main.cc
                 src/fdtd_proxy.cc
                 src/fdtd_solver.cc
                )
  target_include_directories(fdproxy PUBLIC
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                             ${CMAKE_SOURCE_DIR}/external/cxxopts/
                            ${CMAKE_CURRENT_SOURCE_DIR}/../../discretization/fd/fdtd_kernels/include
                            )

  target_link_libraries(fdproxy
    PUBLIC
      discretization
      proxy_utils
      fd_grid
  )

  target_link_kokkos_if_enabled(fdproxy)

endif()
