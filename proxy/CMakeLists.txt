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
