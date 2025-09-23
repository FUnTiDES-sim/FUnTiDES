if(COMPILE_FD)

  add_executable(fdproxy 
                 src/main.cpp
                 src/fdtd_proxy.cpp
                 src/fdtd_solver.cpp
                )
  target_include_directories(fdproxy PUBLIC
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                             ${CMAKE_SOURCE_DIR}/external/cxxopts/
                            ${CMAKE_CURRENT_SOURCE_DIR}/../../discretization/fd/fdtd_kernels/include
                            )

  target_link_libraries(fdproxy
    PRIVATE
      discretization
      proxy_utils
  )

  target_link_kokkos_if_enabled(fdproxy)

endif()