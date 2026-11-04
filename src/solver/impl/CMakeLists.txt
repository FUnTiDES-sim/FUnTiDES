# The fun_solver library: declared here, filled by impl/common (factory and
# shared per-physics types) and by one subtree per compiled method. The
# method selection lives here: COMPILE_SEM/DG/DG_SEM/DG_PADAPTIVE select the
# subtrees and are exported as compile definitions, so the factory and the
# pywrap bindings can guard their per-method parts.
add_library(fun_solver)
add_library(funtides::solver ALIAS fun_solver)
set_target_properties(fun_solver PROPERTIES EXPORT_NAME solver)

target_compile_definitions(fun_solver PUBLIC
    MAX_SOLVER_ACOUSTIC_ORDER=${TARGET_SOLVER_ACOUSTIC_ORDER}
    MAX_SOLVER_ELASTIC_ORDER=${TARGET_SOLVER_ELASTIC_ORDER}
    MAX_SOLVER_ELASTOACOUSTIC_ORDER=${TARGET_SOLVER_ELASTOACOUSTIC_ORDER}
)

target_include_directories(fun_solver PUBLIC
    $<INSTALL_INTERFACE:include/funtides/solver>
)

if(COMPILE_SEM)
  target_compile_definitions(fun_solver PUBLIC COMPILE_SEM)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

# Shared across methods: construction entry point and per-physics types.
add_subdirectory(common)

if(COMPILE_SEM)
  add_subdirectory(sem)
endif()
if(COMPILE_DG)
  add_subdirectory(dg)
endif()
if(COMPILE_DG_SEM)
  add_subdirectory(dgsem)
endif()
if(COMPILE_DG_PADAPTIVE)
  add_subdirectory(dg_padaptive)
endif()

set(SOLV_LIBS fun_discretization fun_core fun_parallel fun_model fun_solver_api)

target_link_kokkos_if_enabled(fun_solver)

target_link_libraries(fun_solver PUBLIC ${SOLV_LIBS})
if(USE_MPI)
target_link_libraries(fun_solver INTERFACE MPI::MPI_CXX)
endif()
