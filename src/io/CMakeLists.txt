add_library(fun_io_api INTERFACE)
add_library(funtides::io_api ALIAS fun_io_api)
set_target_properties(fun_io_api PROPERTIES EXPORT_NAME io_api)

target_include_directories(fun_io_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api/include>
    $<INSTALL_INTERFACE:include/funtides/io>
)

target_link_libraries(fun_io_api
  INTERFACE Kokkos::kokkos fun_core fun_parallel
)

install_funtides_headers(io ${CMAKE_CURRENT_SOURCE_DIR}/api/include)

add_subdirectory(impl)
