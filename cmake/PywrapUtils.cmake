#-------------------------------------------------------------------
# Python wrapping related functions
#-------------------------------------------------------------------

function(configure_python_module target_name)
  set_target_properties(${target_name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

  # set rpath to find proxys shared libraries at runtime
  set_target_properties(${target_name} PROPERTIES
    INSTALL_RPATH "$ORIGIN/../lib64;$ORIGIN/../lib"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
endfunction()

function(install_pyproxys_package)
  # Set install directory for pyproxys
  set(PYPROXYS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/pyproxys")

  # Create the directory at install time
  install(CODE "file(MAKE_DIRECTORY \"${PYPROXYS_INSTALL_DIR}\")")

  # Install libraries from input targets into the directory
  install(TARGETS ${ARGV}
    LIBRARY DESTINATION pyproxys
    ARCHIVE DESTINATION pyproxys
    RUNTIME DESTINATION pyproxys
  )

  # Create __init__.py to make it a Python package
  install(CODE "file(WRITE \"${PYPROXYS_INSTALL_DIR}/__init__.py\" \"\")")
endfunction()