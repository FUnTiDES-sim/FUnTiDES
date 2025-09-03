add_library(model_lib
    src/model_unstruct.cc
    src/model_struct.cc
)

target_include_directories(model_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

set(MESH_LIBS utils)

if (USE_KOKKOS)
    list(APPEND MESH_LIBS Kokkos::kokkos)
endif()

target_link_libraries(model_lib PUBLIC ${MESH_LIBS})
