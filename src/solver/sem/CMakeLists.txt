# SEM solver tree. The fun_solver library target is declared here — every
# variant (dg/, dgsem/, dg_padaptive/) compiles into it — with the
# physics-generic skeleton in common/ and one subdirectory per physics,
# holding that physics' types (wavefield/rhs/traits, impl fragment) beside
# its explicit instantiations.
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

include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

add_subdirectory(common)
add_subdirectory(acoustic)
add_subdirectory(elastic)
add_subdirectory(acoustoelastic)
add_subdirectory(pywrap)

set(SOLV_LIBS fun_discretization fun_core fun_parallel fun_model fun_solver_api)

target_link_kokkos_if_enabled(fun_solver)

target_link_libraries(fun_solver PUBLIC ${SOLV_LIBS})
if(USE_MPI)
target_link_libraries(fun_solver INTERFACE MPI::MPI_CXX)
endif()
