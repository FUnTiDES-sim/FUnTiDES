# Crée la bibliothèque solver
add_library(solver
  src/SEMsolver.cpp
)

# Inclut les headers publics
target_include_directories(solver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Lie avec la bibliothèque utils uniquement
target_link_libraries(solver
  PUBLIC utils mesh Kokkos::kokkos
)
