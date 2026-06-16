target_sources(fun_model PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian_struct_builder.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian_unstruct_builder.cc
)

target_include_directories(fun_model PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
