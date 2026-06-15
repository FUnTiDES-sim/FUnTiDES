target_sources(fun_model PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/model_unstruct.cc
)

target_include_directories(fun_model PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
