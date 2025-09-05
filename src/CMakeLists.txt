add_subdirectory(utils)
add_subdirectory(discretization)
add_subdirectory(solver)
add_subdirectory(model)
add_subdirectory(main)

if(ENABLE_PYWRAP)
  install_pyproxys_package(model solver)
endif()