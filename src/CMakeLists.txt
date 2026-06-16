add_subdirectory(utils)
add_subdirectory(discretization)
add_subdirectory(model)
add_subdirectory(solver)
add_subdirectory(gradient)
add_subdirectory(main)

add_library(funtides INTERFACE)
add_library(funtides::funtides ALIAS funtides)

target_link_libraries(funtides INTERFACE
  fun_utils
  fun_disctretization
  fun_model
  fun_solver
  fun_gradient
)

install(TARGETS funtides fun_utils fun_discretization fun_model fun_solver fun_gradient
        EXPORT funtidesTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
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
