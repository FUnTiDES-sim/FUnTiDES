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
  message(STATUS "")
  
  message(STATUS "Programming Models:")
  message(STATUS "  USE_VECTOR:           ${USE_VECTOR}")
  message(STATUS "  USE_KOKKOS:           ${USE_KOKKOS}")
  message(STATUS "  ENABLE_CUDA:          ${ENABLE_CUDA}")
  message(STATUS "  USE_KOKKOS_TEAMS:     ${USE_KOKKOS_TEAMS}")
  message(STATUS "")

  message(STATUS "Python Wrapping:")
  message(STATUS "  ENABLE_PYWRAP:        ${ENABLE_PYWRAP}")
  message(STATUS "")
  
  message(STATUS "Debugging Options:")
  message(STATUS "  FD_SAVE_SNAPSHOTS:    ${FD_SAVE_SNAPSHOTS}")
  message(STATUS "  PRINT_ALLOC_INFO:     ${PRINT_ALLOC_INFO}")
  message(STATUS "")
  
  message(STATUS "Build Options:")
  message(STATUS "  CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}")
  message(STATUS "  BUILD_SHARED_LIBS:    ${BUILD_SHARED_LIBS}")
  message(STATUS "  CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
  message(STATUS "")
  message(STATUS "==========================================")
  message(STATUS "")
endfunction()