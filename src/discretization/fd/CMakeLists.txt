# FDTD Discretization - Stencils, kernels, and boundary conditions

add_library(fd_discretization INTERFACE)

target_include_directories(fd_discretization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/stencils/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/kernels/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/abckernels/include>
    $<INSTALL_INTERFACE:include>
)

# Dependencies - also use INTERFACE
target_include_directories(fd_discretization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/utils/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/utils/include>
)

# Installation rules
install(TARGETS fd_discretization
    EXPORT fd_discretization-targets
)

install(DIRECTORY stencils/include/
                  kernels/include/
                  boundary/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)
