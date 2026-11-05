#-------------------------------------------------------------------
# Environment script: <install prefix>/funtides_env.sh
#-------------------------------------------------------------------
# Generates a script to source (e.g. from ~/.bashrc) that sets up the TPL
# environment and makes the installed FUnTiDES visible.

# The TPL install prefix is found from Kokkos: <tpl>/lib/cmake/Kokkos
if(NOT DEFINED FUNTIDES_TPL_ENV_SCRIPT AND Kokkos_DIR)
  get_filename_component(_tpl_prefix "${Kokkos_DIR}/../../.." ABSOLUTE)
  if(EXISTS "${_tpl_prefix}/setup_env.sh")
    set(_tpl_env "${_tpl_prefix}/setup_env.sh")
  else()
    set(_tpl_env "")
  endif()
  set(FUNTIDES_TPL_ENV_SCRIPT "${_tpl_env}" CACHE FILEPATH
      "TPL environment script sourced by funtides_env.sh (empty: none)")
endif()

set(FUNTIDES_ENV_PYTHONPATH "")
if(ENABLE_PYWRAP)
  # Keep the python dir relative to the prefix when possible (relocatable)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_PYTHONDIR}")
    file(RELATIVE_PATH _pydir "${CMAKE_INSTALL_PREFIX}" "${CMAKE_INSTALL_PYTHONDIR}")
    if(_pydir MATCHES "^\\.\\.")
      set(_pydir "${CMAKE_INSTALL_PYTHONDIR}")
    else()
      set(_pydir "\${FUNTIDES_INSTALL_DIR}/${_pydir}")
    endif()
  else()
    set(_pydir "\${FUNTIDES_INSTALL_DIR}/${CMAKE_INSTALL_PYTHONDIR}")
  endif()
  set(FUNTIDES_ENV_PYTHONPATH "_funtides_prepend PYTHONPATH \"${_pydir}\"")
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/funtides_env.sh.in
               ${CMAKE_BINARY_DIR}/funtides_env.sh @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/funtides_env.sh DESTINATION .)
