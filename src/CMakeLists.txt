add_subdirectory(utils)
add_subdirectory(discretization)
add_subdirectory(model)
add_subdirectory(solver)
add_subdirectory(main)

if(ENABLE_PYWRAP)
  install_pyfuntides_package(model solver)
endif()


# Collect all sources
file(GLOB_RECURSE ALL_SOURCE_FILES
     "${CMAKE_SOURCE_DIR}/src/*.cpp"
     "${CMAKE_SOURCE_DIR}/src/*.hpp"
     "${CMAKE_SOURCE_DIR}/src/*.cxx"
     "${CMAKE_SOURCE_DIR}/src/*.hxx"
     "${CMAKE_SOURCE_DIR}/src/*.cc"
     "${CMAKE_SOURCE_DIR}/src/*.c"
     "${CMAKE_SOURCE_DIR}/src/*.h"
     "${CMAKE_SOURCE_DIR}/src/*.hh")

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
