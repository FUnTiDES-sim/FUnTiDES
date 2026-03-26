add_library(proxy_solver)

add_subdirectory(acoustic)
add_subdirectory(acoustoelastic)
add_subdirectory(common)
add_subdirectory(elastic)

set(SOLV_LIBS discretization proxy_utils proxy_model_struct proxy_model_unstruct proxy_solver_api proxy_model_api)

target_link_kokkos_if_enabled(proxy_solver)

target_link_libraries(proxy_solver PUBLIC ${SOLV_LIBS})
if(USE_MPI)
target_link_libraries(proxy_solver INTERFACE MPI::MPI_CXX)
endif()

install(TARGETS proxy_solver
  EXPORT proxy_solver_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
