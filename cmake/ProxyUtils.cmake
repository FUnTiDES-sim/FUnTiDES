#-------------------------------------------------------------------
# Proxy app related functions
#-------------------------------------------------------------------

function(print_configuration_summary)
  message(STATUS "")
  message(STATUS "====== PROXY CONFIGURATION SUMMARY ======")
  message(STATUS "")

  message(STATUS "Discretization Methods:")
  message(STATUS "  COMPILE_SEM:          ${COMPILE_SEM}")
  message(STATUS "  COMPILE_FD:           ${COMPILE_FD}")
  message(STATUS "  COMPILE_DG:           ${COMPILE_DG}")
  message(STATUS "")

  message(STATUS "Python Wrapping:")
  message(STATUS "  ENABLE_PYWRAP:        ${ENABLE_PYWRAP}")
  message(STATUS "")

  message(STATUS "Coverage")
  message(STATUS "  ENABLE_COVERAGE:      ${ENABLE_COVERAGE}")
  message(STATUS "")

  message(STATUS "MPI Options:")
  message(STATUS "  USE_MPI:              ${USE_MPI}")
  message(STATUS "")

  message(STATUS "Order generation:")
  message(STATUS "  SEM Solvers:          ${MAX_SOLVER_ORDER}")
  message(STATUS "  Gradient:             ${MAX_DIFFERENTIATOR_ORDER}")
  message(STATUS "")

  message(STATUS "Build Options:")
  message(STATUS "  CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}")
  message(STATUS "  BUILD_SHARED_LIBS:    ${BUILD_SHARED_LIBS}")
  message(STATUS "  CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
  message(STATUS "")
  message(STATUS "==========================================")
  message(STATUS "")
endfunction()
