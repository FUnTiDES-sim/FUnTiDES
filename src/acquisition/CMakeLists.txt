add_library(fun_acquisition INTERFACE)
add_library(funtides::acquisition ALIAS fun_acquisition)
set_target_properties(fun_acquisition PROPERTIES EXPORT_NAME acquisition)

target_include_directories(fun_acquisition
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/funtides/acquisition>
)

install_funtides_headers(acquisition ${CMAKE_CURRENT_SOURCE_DIR}/include)

# source_and_receiver_utils.h evaluates basis functions at source and receiver
# positions, hence the dependency on fun_discretization.
target_link_libraries(fun_acquisition
  INTERFACE
    fun_core
    fun_discretization
)
