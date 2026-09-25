add_library(fun_core INTERFACE)
add_library(funtides::core ALIAS fun_core)
set_target_properties(fun_core PROPERTIES EXPORT_NAME core)

target_include_directories(fun_core
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/core/include>
    $<INSTALL_INTERFACE:include/funtides/core>
)

# data_type_kokkos.h defines the array types on top of Kokkos views.
target_link_libraries(fun_core INTERFACE Kokkos::kokkos)

install_funtides_headers(core ${CMAKE_CURRENT_SOURCE_DIR}/include)
install_funtides_headers(core ${CMAKE_BINARY_DIR}/src/core/include)
