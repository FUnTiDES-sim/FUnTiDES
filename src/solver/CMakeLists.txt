add_library( solver 
             src/solverFactory.cpp
             src/SEMSolver_20.cpp
             src/SEMSolver_21.cpp
             src/SEMSolver_22.cpp
             src/SEMSolver_23.cpp
            )

target_include_directories(
  solver
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         ${CMAKE_CURRENT_SOURCE_DIR}/../discretization
         ${CMAKE_CURRENT_SOURCE_DIR}/../model)

set(SOLV_LIBS discretization mesh utils)

if(USE_KOKKOS)
  list(APPEND SOLV_LIBS Kokkos::kokkos)
endif()

if(USE_OMP)
  find_package(OpenMP REQUIRED)
  list(APPEND SOLV_LIBS OpenMP::OpenMP_CXX)
endif()

target_link_libraries(solver PUBLIC ${SOLV_LIBS})

if(ENABLE_PYWRAP)
  pybind11_add_module(pysolver src/pysolver.cpp)

  target_include_directories(pysolver PRIVATE include)

  target_link_libraries(pysolver PUBLIC ${SOLV_LIBS})

  set_target_properties(pysolver PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
  set_target_properties(pysolver PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()
