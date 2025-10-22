#-------------------------------------------------------------------
# Python wrapping related functions
#-------------------------------------------------------------------

function(configure_python_module target_name)
  set_target_properties(${target_name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

  # set rpath to find funtides shared libraries at runtime
  set_target_properties(${target_name} PROPERTIES
    INSTALL_RPATH "$ORIGIN/../../lib64;$ORIGIN/../../lib"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
endfunction()

function(install_pyfuntides_package)
  # Set install directory for pyfuntides
  set(PYFUNTIDES_INSTALL_DIR "${CMAKE_INSTALL_PYTHONDIR}/pyfuntides")

  # Create the directory at install time
  install(CODE "file(MAKE_DIRECTORY \"${PYFUNTIDES_INSTALL_DIR}\")")

  # Install libraries from input targets into the directory
  install(TARGETS ${ARGV}
    LIBRARY DESTINATION ${PYFUNTIDES_INSTALL_DIR}
    ARCHIVE DESTINATION ${PYFUNTIDES_INSTALL_DIR}
    RUNTIME DESTINATION ${PYFUNTIDES_INSTALL_DIR}
  )

  # Create __init__.py to make it a Python package
  install(CODE "file(WRITE \"${PYFUNTIDES_INSTALL_DIR}/__init__.py\" \"\")")
endfunction()