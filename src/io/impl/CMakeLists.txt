add_library(fun_io
  src/posix_io_controller.cc
  src/io_controller_factory.cc
)
add_library(funtides::io ALIAS fun_io)

target_include_directories(fun_io
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(fun_io
  PUBLIC fun_io_api
)
