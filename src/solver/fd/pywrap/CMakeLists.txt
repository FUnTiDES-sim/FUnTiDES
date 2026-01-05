if(ENABLE_PYWRAP)
  pybind11_add_module(fd_solver_py src/bindings.cpp)

  target_include_directories(fd_solver_py
    PRIVATE
      include
      ${CMAKE_SOURCE_DIR}/src/solver/fd/include
      ${CMAKE_SOURCE_DIR}/src/main/fd/include
      ${CMAKE_SOURCE_DIR}/src/discretization/fd/kernels/include
      ${CMAKE_SOURCE_DIR}/src/discretization/fd/stencils/include
      ${CMAKE_SOURCE_DIR}/src/discretization/fd/abckernels/include
      ${CMAKE_SOURCE_DIR}/src/model/grid/include
      ${CMAKE_SOURCE_DIR}/src/utils/include
      ${CMAKE_BINARY_DIR}/src/utils/include
      ${CMAKE_SOURCE_DIR}/external/cxxopts
      ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include
  )

  target_link_libraries(fd_solver_py
    PUBLIC
      fd_solver
      fd_grid
      proxy_utils
  )

  configure_python_module(fd_solver_py)

  # Set output name to fd_solver (without _py suffix)
  set_target_properties(fd_solver_py PROPERTIES OUTPUT_NAME "fd_solver")
endif()
