add_library(fun_acquisition INTERFACE)
add_library(funtides::acquisition ALIAS fun_acquisition)

target_include_directories(fun_acquisition
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# source_and_receiver_utils.h evaluates basis functions at source and receiver
# positions, hence the dependency on fun_discretization.
target_link_libraries(fun_acquisition
  INTERFACE
    fun_core
    fun_discretization
)
