add_subdirectory(src/makutu)

if (ENABLE_Shiva)
  add_subdirectory(src/shiva)
endif()

target_include_directories(proxy_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)