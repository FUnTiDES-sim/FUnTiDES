add_library(solver
  src/SEMsolver.cpp
)

target_include_directories(solver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

set(SOLV_LIBS mesh utils)

if(USE_KOKKOS)
  list(APPEND SOLV_LIBS Kokkos::kokkos)
endif()

if(USE_OMP)
  # find_package(OpenMP REQUIRED)
  list(APPEND SOLV_LIBS OpenMP::OpenMP_CXX)
endif()

target_link_libraries(solver
  PUBLIC ${SOLV_LIBS}
)
