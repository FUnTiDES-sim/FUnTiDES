add_subdirectory(impl)

if(COMPILE_DG)
add_subdirectory(dg)
endif()

if(COMPILE_DG_SEM)
add_subdirectory(dg-sem)
endif()

if(COMPILE_DG_PADAPTIVE)
add_subdirectory(dg_padaptive)
endif()
