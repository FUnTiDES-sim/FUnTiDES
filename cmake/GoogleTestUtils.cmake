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

    gtest_discover_tests(${TEST_NAME})
endfunction()