add_subdirectory(utils)
add_subdirectory(discretization)
add_subdirectory(model)
add_subdirectory(solver)
add_subdirectory(main)

if(ENABLE_PYWRAP)
  install_pyproxys_package(model solver)
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

# Example: filter out files whose path contains "thirdparty" or "external"
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "thirdparty")
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "external")

# Make the format target
add_custom_target( apply_clang_format
  COMMAND clang-format -i
          -style=file:${CMAKE_SOURCE_DIR}/.clang-format
          ${ALL_SOURCE_FILES}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
