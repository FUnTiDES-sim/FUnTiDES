add_library(proxy_utils INTERFACE)

target_include_directories(proxy_utils
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/utils/include>
)

add_library(proxy_adios INTERFACE)

target_include_directories(proxy_adios INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(proxy_adios INTERFACE adios2::adios2)

target_sources(proxy_adios INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/adios_utils.h>
)
