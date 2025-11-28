add_subdirectory(utils)
add_subdirectory(discretization)
add_subdirectory(model)
add_subdirectory(solver)
add_subdirectory(main)

if(ENABLE_PYWRAP)
  install_pyfuntides_package(model solver)
endif()


# Collect all sources
set(DIRS "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/test")
set(EXTS cpp hpp cxx hxx cc c h hh)

set(GLOBS "")
foreach(d ${DIRS})
  foreach(e ${EXTS})
    message( STATUS "Adding source files from dir: ${d} with extension: .${e}" )
    list(APPEND GLOBS "${d}/*.${e}")
  endforeach()
endforeach()

file(GLOB_RECURSE ALL_SOURCE_FILES ${GLOBS})

# Example: filter out files whose path contains "Shiva" or "SEMKernels"
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "Shiva")
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
