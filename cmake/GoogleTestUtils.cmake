#-------------------------------------------------------------------
# Google test related functions
#-------------------------------------------------------------------

function(add_gtest TEST_NAME TEST_FILE EXTRA_LINKS)
    add_executable(${TEST_NAME}
        ${TEST_FILE}
        ${TEST_MAIN}
    )

    target_link_libraries(${TEST_NAME}
        PRIVATE
        gtest
        gtest_main
        ${EXTRA_LINKS}
    )

    target_link_kokkos_if_enabled(${TEST_NAME})

    add_dependencies(build_tests ${TEST_NAME})

    gtest_discover_tests(${TEST_NAME} DISCOVERY_TIMEOUT 60)
endfunction()

# Group several test .cc files into a single executable: one main.cc
# compile + one link per group instead of one per file. Each source stays
# its own translation unit; only the executable is shared.
function(add_gtest_group GROUP_NAME)
    set(multiValueArgs SOURCES EXTRA_LINKS)
    cmake_parse_arguments(GRP "" "" "${multiValueArgs}" ${ARGN})

    add_executable(${GROUP_NAME}
        ${GRP_SOURCES}
        ${TEST_MAIN}
    )

    target_link_libraries(${GROUP_NAME}
        PRIVATE
        gtest
        gtest_main
        ${GRP_EXTRA_LINKS}
    )

    target_link_kokkos_if_enabled(${GROUP_NAME})

    add_dependencies(build_tests ${GROUP_NAME})

    gtest_discover_tests(${GROUP_NAME} DISCOVERY_TIMEOUT 60)
endfunction()
