add_gtest(semkernels_qk_hexahedron_test
  test_Qk_Hexahedron.cpp
  discretization
)

# Add include directories
target_include_directories(semkernels_qk_hexahedron_test PRIVATE
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/finiteElement
  ${CMAKE_SOURCE_DIR}/src/finiteElement/makutu
  ${CMAKE_SOURCE_DIR}/src/common
)