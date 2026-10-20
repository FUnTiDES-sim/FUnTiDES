# Common boundary of every solver method: the pure interfaces in api/ (base
# classes consumed by SEM, DG, DG-SEM and the p-adaptive variants) and the
# construction entry point in factory/. The method trees (sem/, dg/, dgsem/,
# dg_padaptive/) hang beside them and all compile into fun_solver.
add_subdirectory(api)

if(COMPILE_SEM)
  add_subdirectory(sem)
endif()

# fun_solver is declared in sem/, so the factory (which compiles into it and
# guards DG variants with the COMPILE_DG* definitions) comes after.
add_subdirectory(factory)

if(COMPILE_DG)
  add_subdirectory(dg)
endif()
if(COMPILE_DG_SEM)
  add_subdirectory(dgsem)
endif()
if(COMPILE_DG_PADAPTIVE)
  add_subdirectory(dg_padaptive)
endif()
