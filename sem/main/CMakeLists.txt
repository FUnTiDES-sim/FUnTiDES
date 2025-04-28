add_executable(semproxy
  src/main.cpp
  src/SEMproxy.cpp
)

target_include_directories(semproxy
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(semproxy
  PRIVATE
    solver
    mesh
    utils
)

if (ENABLE_PYWRAP)
  file(COPY ${CMAKE_SOURCE_DIR}/scripts/pyscript/ DESTINATION ${CMAKE_BINARY_DIR}/proxy)
  pybind11_add_module(pysem src/pysem.cpp src/SEMproxy.cpp)

  target_include_directories(pysem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

  target_link_libraries(pysem
    PRIVATE
      solver
      mesh
      utils
  )

  set_target_properties(pysem PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(pysem PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()
