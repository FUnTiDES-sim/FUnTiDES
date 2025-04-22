add_executable(semproxy
  src/main.cpp
  src/SEMproxy.cpp
)

target_include_directories(semproxy
  PRIVATE
    ${CMAKE_SOURCE_DIR}/proxy/include
)

target_link_libraries(semproxy
  PRIVATE
    solver
    mesh
    utils
)

file(COPY ${CMAKE_SOURCE_DIR}/pyscript/ DESTINATION ${CMAKE_BINARY_DIR}/proxy)
pybind11_add_module(pysem src/pysem.cpp src/SEMproxy.cpp)

target_include_directories(pysem PRIVATE ${CMAKE_SOURCE_DIR}/proxy/include)

target_link_libraries(pysem
  PRIVATE
    solver
    mesh
    utils
)
