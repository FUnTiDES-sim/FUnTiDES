#-------------------------------------------------------------------
# Python wrapping related functions
#-------------------------------------------------------------------

# fix warning yield when install path is ./python
if(POLICY CMP0177)
    cmake_policy(SET CMP0177 NEW)
endif()

function(configure_python_module target_name)
  set_target_properties(${target_name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

  # set rpath to find funtides shared libraries at runtime
  if(SKBUILD)
    set(rpath "$ORIGIN")
  else()
    set(rpath "$ORIGIN/../../lib64;$ORIGIN/../../lib")
  endif()
  set_target_properties(${target_name} PROPERTIES
    INSTALL_RPATH "${rpath}"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
endfunction()

# Set install directory for pyfuntides
function(install_pyfuntides_package)
  if(SKBUILD)
    set(PYFUNTIDES_INSTALL_DIR "pyfuntides")
  else()
    set(PYFUNTIDES_INSTALL_DIR "${CMAKE_INSTALL_PYTHONDIR}/pyfuntides")
  endif()

  # Configure _build_info.py with variables only known at build time, such
  # as the compiled methods and max orders
  configure_file( "${PROJECT_SOURCE_DIR}/src/pyfuntides/_build_info.py.in"
                  "${PROJECT_BINARY_DIR}/src/pyfuntides/_build_info.py" @ONLY)

  # Install libraries from input targets into the directory
  install(TARGETS ${ARGV}
    LIBRARY DESTINATION ${PYFUNTIDES_INSTALL_DIR} COMPONENT python
    ARCHIVE DESTINATION ${PYFUNTIDES_INSTALL_DIR} COMPONENT python
    RUNTIME DESTINATION ${PYFUNTIDES_INSTALL_DIR} COMPONENT python
  )

  # Install the package's Python files (__init__.py and the generated _build_info.py)
  install(FILES "${PROJECT_SOURCE_DIR}/src/pyfuntides/__init__.py"
                "${PROJECT_BINARY_DIR}/src/pyfuntides/_build_info.py"
          DESTINATION "${PYFUNTIDES_INSTALL_DIR}"
          COMPONENT python)
endfunction()
