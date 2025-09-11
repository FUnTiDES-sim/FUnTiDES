<<<<<<< HEAD
add_library(solver
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

# Allow to link agains any model and mesh library that respects the API
if(NOT DEFINED SOLVER_MODEL_INTERFACE)
  set(SOLVER_MODEL_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../model)
endif()

if(NOT DEFINED SOLVER_MESH_LIBS)
  set(SOLVER_MESH_LIBS model_lib)
endif()

message(STATUS "SOLVER_MODEL_INTERFACE: ${SOLVER_MODEL_INTERFACE}")
message(STATUS "SOLVER_MESH_LIBS      : ${SOLVER_MESH_LIBS}")

target_include_directories(
  solver
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         ${CMAKE_CURRENT_SOURCE_DIR}/../discretization
         ${SOLVER_MODEL_INTERFACE})

set(SOLV_LIBS discretization ${SOLVER_MESH_LIBS} utils)

if(USE_KOKKOS)
  list(APPEND SOLV_LIBS Kokkos::kokkos)
endif()

if(USE_OMP)
  find_package(OpenMP REQUIRED)
  list(APPEND SOLV_LIBS OpenMP::OpenMP_CXX)
endif()

target_link_libraries(solver PUBLIC ${SOLV_LIBS})

if(ENABLE_PYWRAP)
  pybind11_add_module(pysolver src/pysolver.cpp src/SEMsolver.cpp)

  target_include_directories(pysolver PRIVATE
             include
             ${CMAKE_CURRENT_SOURCE_DIR}/../discretization
             ${SOLVER_MODEL_INTERFACE}
             ${CMAKE_SOURCE_DIR}/external/pykokkos-base/include)

  target_link_libraries(pysolver PUBLIC ${SOLV_LIBS})

  set_target_properties(pysolver PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(pysolver PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()
=======
add_subdirectory(api)
add_subdirectory(impl)
add_subdirectory(pywrap)
>>>>>>> origin/main
