add_library(proxy_solver
            src/solver_factory.cc
            src/geos/SEMQ1_struct_geos_solver.cpp
            src/geos/SEMQ1_unstruct_geos_solver.cpp
            src/geos/SEMQ2_struct_geos_solver.cpp
            src/geos/SEMQ2_unstruct_geos_solver.cpp
            src/geos/SEMQ3_struct_geos_solver.cpp
            src/geos/SEMQ3_unstruct_geos_solver.cpp)
            # src/shiva/SEMQ1_struct_shiva_solver.cpp
            # src/shiva/SEMQ1_unstruct_shiva_solver.cpp
            # src/shiva/SEMQ2_struct_shiva_solver.cpp
            # src/shiva/SEMQ2_unstruct_shiva_solver.cpp
            # src/shiva/SEMQ3_struct_shiva_solver.cpp
            # src/shiva/SEMQ3_unstruct_shiva_solver.cpp)

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
