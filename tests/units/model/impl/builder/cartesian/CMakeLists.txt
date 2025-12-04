add_gtest(cartesian_struct_builder_test
  test_cartesian_struct_builder.cc
  proxy_model_builder_cartesian
)

add_gtest(sep_params_test
  test_sep_params.cc
  proxy_model_builder_cartesian
)

# Copy SEP test data to build directory (relative to where tests run)
file(COPY ${PROJECT_SOURCE_DIR}/tests/units/data/sep/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/units/data/sep/)
