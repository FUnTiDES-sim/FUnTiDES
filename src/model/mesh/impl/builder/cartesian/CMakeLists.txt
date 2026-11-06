target_sources(fun_model PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian_struct_builder.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian_unstruct_builder.cc
)

if(COMPILE_SEM_2D)
  target_compile_definitions(fun_model PUBLIC COMPILE_SEM_2D)
  target_sources(fun_model PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian_unstruct_2d_builder.cc
  )
endif()

target_include_directories(fun_model PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(model ${CMAKE_CURRENT_SOURCE_DIR}/include)
