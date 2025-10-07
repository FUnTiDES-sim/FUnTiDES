add_executable(cartesian_struct_builder_test
    test_cartesian_struct_builder.cc
)

target_link_libraries(cartesian_struct_builder_test
    PRIVATE
    gtest
    gtest_main
    proxy_model_builder_cartesian
)


gtest_discover_tests(cartesian_struct_builder_test)