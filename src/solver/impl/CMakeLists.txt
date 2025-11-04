add_library(proxy_solver)

add_subdirectory(src/makutu)
target_sources(proxy_solver PRIVATE ${MAKUTU_GENERATED_SOURCES})
#add_subdirectory(src/shiva)
#target_sources(proxy_solver PRIVATE ${SHIVA_GENERATED_SOURCES})

target_sources(proxy_solver PRIVATE 
     src/solver_factory.cc
)

target_include_directories(proxy_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

set(SOLV_LIBS discretization proxy_utils proxy_model_struct proxy_model_unstruct proxy_solver_api proxy_model_api)

target_link_kokkos_if_enabled(proxy_solver)

target_link_libraries(proxy_solver PUBLIC ${SOLV_LIBS})

install(TARGETS proxy_solver
  EXPORT proxy_solver_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
