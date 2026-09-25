add_library(fun_io_api INTERFACE)
add_library(funtides::io_api ALIAS fun_io_api)

target_include_directories(fun_io_api
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api/include>
)

target_link_libraries(fun_io_api
  INTERFACE Kokkos::kokkos fun_core fun_parallel
)

add_subdirectory(impl)
