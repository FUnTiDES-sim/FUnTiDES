if(COMPILE_FD)

  add_executable(fdproxy src/FDTDmain.cpp)

  target_include_directories(fdproxy
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/../discretization
  )

  target_link_libraries(fdproxy
    PRIVATE
      discretization
      proxy_utils
  )

  target_link_kokkos_if_enabled(fdproxy)

endif()