add_library(cartesian_builders
    src/cartesian_struct_mesh.cc
    src/cartesian_unstruct_mesh.cc
)

target_include_directories(cartesian_builders PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../model_api/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../model_builder_api/include>
)

set(BUILDER_LIBS utils model_lib)

if (USE_KOKKOS)
    list(APPEND BUILDER_LIBS Kokkos::kokkos)
endif()

target_link_libraries(cartesian_builders PUBLIC ${BUILDER_LIBS})
