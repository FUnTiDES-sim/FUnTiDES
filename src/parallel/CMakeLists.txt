add_library(fun_parallel INTERFACE)
add_library(funtides::parallel ALIAS fun_parallel)

target_include_directories(fun_parallel
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# topology_factory.h inspects a ModelApi to discover partition boundaries, hence
# the dependency on fun_model. The reverse edge does not exist: no business
# component knows about partitioning.
target_link_libraries(fun_parallel
  INTERFACE
    fun_core
    fun_model
)
