add_library(funtides_io
  src/posix_io_controller.cc
  src/io_controller_factory.cc
)

target_include_directories(funtides_io
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(funtides_io
  PUBLIC funtides_io_api
)
