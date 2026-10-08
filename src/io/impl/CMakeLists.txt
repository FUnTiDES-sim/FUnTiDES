add_library(fun_io
  src/posix_io_controller.cc
  src/io_controller_factory.cc
)
add_library(funtides::io ALIAS fun_io)
set_target_properties(fun_io PROPERTIES EXPORT_NAME io)

target_include_directories(fun_io
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(fun_io
  PUBLIC fun_io_api
)

if(FUNTIDES_ENABLE_ADIOS2)
  # 2.9 is the first release with adios2::Mode::ReadRandomAccess.
  find_package(ADIOS2 2.9 REQUIRED)
  find_package(Threads REQUIRED)

  # ADIOS2 >= 2.11 renamed the C++ bindings target to adios2::cxx and marks
  # adios2::cxx11 as deprecated; older releases only provide adios2::cxx11.
  if(TARGET adios2::cxx)
    set(_fun_adios2_cxx adios2::cxx)
  else()
    set(_fun_adios2_cxx adios2::cxx11)
  endif()

  target_sources(fun_io PRIVATE src/adios2_io_controller.cc)
  # Serial C++ bindings on purpose: each rank writes its own BP file and the
  # I/O thread must never issue MPI calls (see adios2_io_controller.h).
  # PRIVATE: no ADIOS2 type appears in the public API.
  target_link_libraries(fun_io PRIVATE ${_fun_adios2_cxx} Threads::Threads)
  target_compile_definitions(fun_io PRIVATE FUNTIDES_IO_HAVE_ADIOS2)
  message(STATUS "FUnTiDES io: ADIOS2 backend enabled (ADIOS2 ${ADIOS2_VERSION})")
else()
  message(STATUS "FUnTiDES io: ADIOS2 backend disabled")
endif()
