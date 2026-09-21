add_library(funtides_io_api INTERFACE)

target_include_directories(funtides_io_api
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/api/include
)

target_link_libraries(funtides_io_api
  INTERFACE Kokkos::kokkos fun_utils
)

add_subdirectory(impl)
