add_subdirectory(core)
add_subdirectory(discretization)
add_subdirectory(model)
add_subdirectory(parallel)
add_subdirectory(acquisition)
add_subdirectory(solver)
add_subdirectory(gradient)
add_subdirectory(io)
add_subdirectory(main)

add_library(funtides INTERFACE)
add_library(funtides::funtides ALIAS funtides)

target_link_libraries(funtides INTERFACE
  fun_core
  fun_parallel
  fun_acquisition
  fun_discretization
  fun_model
  fun_solver
  fun_gradient
  fun_io
)

install(TARGETS funtides
                fun_core
                fun_parallel
                fun_acquisition
                fun_discretization fun_discretization_api fun_discretization_impl
                fun_discretization_common fun_discretization_makutu
                fun_discretization_tensorial
                fun_model
                fun_solver_api fun_solver
                fun_gradient_api fun_gradient
                fun_io_api fun_io
        EXPORT funtidesTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
)

# CMake package: find_package(funtides) then link funtides::funtides, or a
# single component such as funtides::solver.
include(CMakePackageConfigHelpers)

install(EXPORT funtidesTargets
        NAMESPACE funtides::
        DESTINATION lib/cmake/funtides
)

configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/funtidesConfig.cmake.in
  ${CMAKE_BINARY_DIR}/funtidesConfig.cmake
  INSTALL_DESTINATION lib/cmake/funtides
)

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/funtidesConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES ${CMAKE_BINARY_DIR}/funtidesConfig.cmake
              ${CMAKE_BINARY_DIR}/funtidesConfigVersion.cmake
        DESTINATION lib/cmake/funtides
)

if(ENABLE_PYWRAP)
  install_pyfuntides_package(model solver gradient)
endif()


## Clang formating
# Collect all sources
set(DIRS "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/tests")
set(EXTS cpp hpp cxx hxx cc c h hh)

set(GLOBS "")
foreach(d ${DIRS})
  foreach(e ${EXTS})
    message( STATUS "Adding source files from dir: ${d} with extension: .${e}" )
    list(APPEND GLOBS "${d}/*.${e}")
  endforeach()
endforeach()

file(GLOB_RECURSE ALL_SOURCE_FILES ${GLOBS})

# Example: filter out files whose path contains "SEMKernels"
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "SEMKernels")

# Make the format target
if( NOT CLANG_FORMAT_EXECUTABLE )
  set( CLANG_FORMAT_EXECUTABLE clang-format )
endif()

message( STATUS "Using clang-format executable: ${CLANG_FORMAT_EXECUTABLE}" )

add_custom_target( apply_clang_format
                   COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file -fallback-style=none ${ALL_SOURCE_FILES}
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
