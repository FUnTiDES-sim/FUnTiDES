#-------------------------------------------------------------------
# Python wrapping related functions
#-------------------------------------------------------------------

function(configure_python_module target_name)
  set_target_properties(${target_name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endfunction()