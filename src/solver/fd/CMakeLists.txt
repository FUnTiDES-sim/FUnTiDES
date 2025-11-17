# FD Solver Library
# Time integration and source/receiver management for FDTD simulations

add_library(fd_solver
    src/fd_solver.cc
)

target_include_directories(fd_solver PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Add include paths for dependencies
target_include_directories(fd_solver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/main/fd/include
    ${CMAKE_SOURCE_DIR}/src/utils/include
    ${CMAKE_BINARY_DIR}/src/utils/include
    ${CMAKE_SOURCE_DIR}/src/discretization/fd/kernels/include
    ${CMAKE_SOURCE_DIR}/src/discretization/fd/stencils/include
    ${CMAKE_SOURCE_DIR}/src/discretization/fd/abckernels/include
    ${CMAKE_SOURCE_DIR}/src/model/grid/include
    ${CMAKE_SOURCE_DIR}/external/cxxopts
)

target_link_libraries(fd_solver PUBLIC
    fd_grid
)

# Optional: Add Kokkos support if available
if(USE_KOKKOS)
    target_link_libraries(fd_solver PUBLIC Kokkos::kokkos)
endif()

# Installation rules
install(TARGETS fd_solver
    EXPORT fd_solver-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)
