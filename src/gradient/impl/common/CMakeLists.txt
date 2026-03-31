target_sources(funtides_gradient PRIVATE 
     src/differentiator_factory.cc
)

target_include_directories(funtides_gradient
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
