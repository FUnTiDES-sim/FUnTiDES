add_executable(validate_solution
    validate_solution.cc
)

# ==============================================================================
# Path configuration (ABSOLUTE paths to work with any build directory name)
# ==============================================================================

# These paths will be injected as compile-time definitions
# This ensures the validation works regardless of the build directory name
# Reference files are stored in tests/analyticalsolution/
target_compile_definitions(validate_solution PRIVATE
    FUNTIDES_EXECUTABLE="${CMAKE_BINARY_DIR}/bin/funtidessem"
    EXTRACT_RECEIVERS_SCRIPT="${CMAKE_SOURCE_DIR}/scripts/adios/adios_single_receiver_viz.py"
    VALIDATION_REFERENCE_ACOUSTIC="${CMAKE_SOURCE_DIR}/tests/validation/analyticalsolution/P.dat"
    VALIDATION_REFERENCE_ELASTIC="${CMAKE_SOURCE_DIR}/tests/validation/analyticalsolution/Ux.dat"
)

# ==============================================================================
# Validation targets
# ==============================================================================

# Configurable tolerance variables (separate for acoustic and elastic)
set(VALIDATION_TOLERANCE_ACOUSTIC "5e-2" CACHE STRING "Acoustic validation tolerance (5%)")
set(VALIDATION_TOLERANCE_ELASTIC "0.2" CACHE STRING "Elastic validation tolerance (20%)")

# Main target: make validate
# Runs all validation tests
add_custom_target(validate
    COMMAND ${CMAKE_CTEST_COMMAND} -L validation --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running all validation tests (24 tests: 2 meshes × 3 orders × 2 types × 2 model-on-nodes)"
)

# Acoustic only tests
add_custom_target(validate_acoustic
    COMMAND ${CMAKE_CTEST_COMMAND} -L "validation_acoustic" --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running acoustic validation tests"
)

# Elastic only tests
add_custom_target(validate_elastic
    COMMAND ${CMAKE_CTEST_COMMAND} -L "validation_elastic" --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running elastic validation tests"
)

# Cartesian mesh only tests
add_custom_target(validate_cartesian
    COMMAND ${CMAKE_CTEST_COMMAND} -L "mesh_cartesian" --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running cartesian mesh validation tests"
)

# Ucartesian mesh only tests
add_custom_target(validate_ucartesian
    COMMAND ${CMAKE_CTEST_COMMAND} -L "mesh_ucartesian" --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running ucartesian mesh validation tests"
)

# Model on nodes tests
add_custom_target(validate_model_on_nodes
    COMMAND ${CMAKE_CTEST_COMMAND} -L "model_on_nodes" --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running model-on-nodes validation tests"
)

# Model NOT on nodes tests
add_custom_target(validate_no_model_on_nodes
    COMMAND ${CMAKE_CTEST_COMMAND} -L "no_model_on_nodes" --output-on-failure
    DEPENDS validate_solution funtidessem
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running validation tests without model-on-nodes"
)

# Cleanup target: make clean_validation
add_custom_target(clean_validation
    COMMAND ${CMAKE_COMMAND} -E remove receivers_numerical.txt
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Cleaning validation output files"
)

# ==============================================================================
# CTest integration - Parameterized tests
# ==============================================================================

# Loop over all configuration combinations
foreach(MESH_TYPE cartesian ucartesian)
    foreach(ORDER 1 2 3)
        foreach(MODEL_ON_NODES OFF ON)
            # Set model-on-nodes flag and label
            if(MODEL_ON_NODES)
                set(MODEL_FLAG "--is-model-on-nodes")
                set(MODEL_LABEL "model_on_nodes")
                set(MODEL_SUFFIX "_model_on_nodes")
            else()
                set(MODEL_FLAG "")
                set(MODEL_LABEL "no_model_on_nodes")
                set(MODEL_SUFFIX "")
            endif()
            
            # Acoustic test for this configuration
            add_test(
                NAME validation_acoustic_order${ORDER}_${MESH_TYPE}${MODEL_SUFFIX}
                COMMAND validate_solution 
                    --order ${ORDER}
                    --mesh ${MESH_TYPE}
                    ${MODEL_FLAG}
                    --tolerance ${VALIDATION_TOLERANCE_ACOUSTIC}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            )
            
            set_tests_properties(validation_acoustic_order${ORDER}_${MESH_TYPE}${MODEL_SUFFIX} PROPERTIES
                LABELS "validation;validation_acoustic;order${ORDER};mesh_${MESH_TYPE};${MODEL_LABEL}"
                TIMEOUT 600
            )
            
            # Elastic test for this configuration
            add_test(
                NAME validation_elastic_order${ORDER}_${MESH_TYPE}${MODEL_SUFFIX}
                COMMAND validate_solution 
                    --order ${ORDER}
                    --elastic
                    --mesh ${MESH_TYPE}
                    ${MODEL_FLAG}
                    --tolerance ${VALIDATION_TOLERANCE_ELASTIC}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            )
            
            set_tests_properties(validation_elastic_order${ORDER}_${MESH_TYPE}${MODEL_SUFFIX} PROPERTIES
                LABELS "validation;validation_elastic;order${ORDER};mesh_${MESH_TYPE};${MODEL_LABEL}"
                TIMEOUT 600
            )
        endforeach()
    endforeach()
endforeach()
