# FDTD Discretization - Stencils, kernels, and boundary conditions

add_library(fdtd_discretization INTERFACE)

target_include_directories(fdtd_discretization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/stencils/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/kernels/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/boundary/include>
    $<INSTALL_INTERFACE:include>
)

# Dependencies - also use INTERFACE
target_include_directories(fdtd_discretization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/utils/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/utils/include>
)

# Installation rules
install(TARGETS fdtd_discretization
    EXPORT fdtd_discretization-targets
)

install(DIRECTORY stencils/include/
                  kernels/include/
                  boundary/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)
