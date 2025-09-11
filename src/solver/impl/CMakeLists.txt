add_library(proxy_solver
            src/solverFactory.cpp
            src/SEMQ1_classic_solver.cpp
            src/SEMQ2_classic_solver.cpp
            src/SEMQ3_classic_solver.cpp
            src/SEMQ1_optim_solver.cpp
            src/SEMQ2_optim_solver.cpp
            src/SEMQ3_optim_solver.cpp
            src/SEMQ1_geos_solver.cpp
            src/SEMQ2_geos_solver.cpp
            src/SEMQ3_geos_solver.cpp
            src/SEMQ1_shiva_solver.cpp
            src/SEMQ2_shiva_solver.cpp
            src/SEMQ3_struct_shiva_solver.cpp
            src/SEMQ3_unstruct_shiva_solver.cpp)

target_include_directories(proxy_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

set(SOLV_LIBS discretization proxy_utils proxy_model_struct proxy_model_unstruct proxy_solver_api proxy_model_api)

target_link_kokkos_if_enabled(proxy_solver)

if(USE_OMP)
  find_package(OpenMP REQUIRED)
  list(APPEND SOLV_LIBS OpenMP::OpenMP_CXX)
endif()

target_link_libraries(proxy_solver PUBLIC ${SOLV_LIBS})

install(TARGETS proxy_solver
  EXPORT proxy_solver_targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
