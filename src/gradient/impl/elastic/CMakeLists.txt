if(COMPILE_ELASTIC)
  add_subdirectory(src/makutu)
endif()

target_include_directories(fun_gradient
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(gradient ${CMAKE_CURRENT_SOURCE_DIR}/include)
