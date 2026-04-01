#
# Specify list of tests
#

set( unit_tests_sources
        testSEMQkGLIntegrals.cpp
     )

set( dependencyList ${parallelDeps} gtest finiteElement)

#
# Add gtest C++ based tests
#
foreach(test ${unit_tests_sources})
#   message("test is ${test}")
    set( header ${test} )
    string(REPLACE "test" "../" header ${header})
    string(REPLACE ".cpp" ".hpp" header ${header})
#   message("header is ${header}")

    get_filename_component( test_name ${test} NAME_WE )
    add_executable( ${test_name}
                    ${test} 
                    ${dependencyList}
                   )

    target_include_directories( ${test_name}
                                PUBLIC
                                ${header}
                              )

    target_link_libraries( ${test_name}
                           PRIVATE
                           ${dependencyList}
                          )
endforeach()
