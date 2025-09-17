if(COMPILE_SEM)

  add_executable(semproxy
    src/main.cpp
    src/SEMproxy.cpp
  )

  target_include_directories(semproxy
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${CMAKE_SOURCE_DIR}/external/cxxopts
  )

  target_link_libraries(semproxy
    PRIVATE
      discretization
      proxy_solver
      proxy_model_builder_cartesian
      proxy_model_struct
      proxy_model_unstruct
      proxy_utils
      ${STDCXXFS_LIB}
  )

endif()
