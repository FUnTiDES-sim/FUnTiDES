target_sources(funtides_gradient PRIVATE 
     src/differentiator_factory.cc
)

target_compile_definitions(funtides_gradient PUBLIC
    MAX_DIFFERENTIATOR_ACOUSTIC_ORDER=${TARGET_DIFFERENTIATOR_ACOUSTIC_ORDER}
    MAX_DIFFERENTIATOR_ELASTIC_ORDER=${TARGET_DIFFERENTIATOR_ELASTIC_ORDER}
)

target_include_directories(funtides_gradient
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
