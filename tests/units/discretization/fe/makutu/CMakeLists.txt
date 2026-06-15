# ============================================================================
# Generate one test binary per (polynomial order, test group).
# Each binary compiles only one element type and one subset of tests, enabling
# fully parallel compilation and reducing per-TU template instantiation cost.
# ============================================================================

set(_QK_TEST_GROUPS
  mass_stiffness
  jacobian
  interp_gradient
  indexing_bmatrix
  face_interface
  sumfact_acoustic
  sumfact_elastic
)

set(_QK_TEST_ORDERS "")
foreach(_order RANGE 1 ${TARGET_SOLVER_ACOUSTIC_ORDER})
  list(APPEND _QK_TEST_ORDERS ${_order})
endforeach()

foreach(ORDER ${_QK_TEST_ORDERS})
  set(QK_TYPE "Q${ORDER}_Hexahedron_Lagrange_GaussLobatto")

  foreach(GROUP ${_QK_TEST_GROUPS})
    set(_gen_file "${CMAKE_CURRENT_BINARY_DIR}/test_Qk_Hexahedron_Q${ORDER}_${GROUP}.cpp")

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/templates/test_Qk_Hexahedron_group.cpp.in"
      "${_gen_file}"
      @ONLY
    )
    message(STATUS "Configured: test_Qk_Hexahedron_Q${ORDER}_${GROUP}.cpp")

    set(_target "semkernels_qk_hexahedron_Q${ORDER}_${GROUP}_test")
    add_gtest(${_target}
      ${_gen_file}
      fun_discretization
    )

    target_include_directories(${_target} PRIVATE
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_SOURCE_DIR}/src/finiteElement
      ${CMAKE_SOURCE_DIR}/src/finiteElement/makutu
      ${CMAKE_SOURCE_DIR}/src/common
      ${CMAKE_CURRENT_SOURCE_DIR}   # for impl/ directory
    )
  endforeach()
endforeach()

unset(_QK_TEST_ORDERS)
unset(_QK_TEST_GROUPS)
unset(_gen_file)
unset(_target)

# ============================================================================
# lagrange_basis group runs for ALL orders 1-9 so that LagrangeBasis4GL-9GL
# functions (compiled but never called by existing tests) get coverage.
# ============================================================================

set(GROUP "lagrange_basis")
foreach(ORDER RANGE 1 9)
  set(QK_TYPE "Q${ORDER}_Hexahedron_Lagrange_GaussLobatto")
  set(_gen_file "${CMAKE_CURRENT_BINARY_DIR}/test_Qk_Hexahedron_Q${ORDER}_lagrange_basis.cpp")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/templates/test_Qk_Hexahedron_group.cpp.in"
    "${_gen_file}"
    @ONLY
  )

  set(_target "semkernels_qk_hexahedron_Q${ORDER}_lagrange_basis_test")
  add_gtest(${_target}
    ${_gen_file}
    fun_discretization
  )

  target_include_directories(${_target} PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/finiteElement
    ${CMAKE_SOURCE_DIR}/src/finiteElement/makutu
    ${CMAKE_SOURCE_DIR}/src/common
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
endforeach()

unset(_gen_file)
unset(_target)
