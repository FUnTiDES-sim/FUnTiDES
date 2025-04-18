add_library(mesh
    src/SEMmesh.cpp
)

target_include_directories(mesh PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    # $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/utils/include/utils>
)

target_link_libraries(mesh PUBLIC utils Kokkos::kokkos)
