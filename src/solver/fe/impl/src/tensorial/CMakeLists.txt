# ============================================================================
# Configuration: Define all combinations to generate
# ============================================================================

# Define polynomial orders
set(ORDERS 1 2 3)

# Define mesh types
set(MESH_TYPES "struct" "unstruct")

# Define model types
set(MODEL_TYPES "onelements" "onnodes")

set(PHYSIC_TYPE "acoustic" "elastic")

# ============================================================================
# Setup paths
# ============================================================================

# Output directory for generated files (in build tree)
set(GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}")

# ============================================================================
# Generate solver implementation files
# ============================================================================

set(TENSORIAL_GENERATED_SOURCES)

foreach(ORDER ${ORDERS})
    foreach(MESH_TYPE ${MESH_TYPES})
        foreach(PHYSIC_TYPE ${PHYSIC_TYPE})
            set(TEMPLATE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/templates/tensorial_solver_template_${PHYSIC_TYPE}.cpp.in")
            foreach(MODEL_TYPE ${MODEL_TYPES})
                # Create filename: SEMQ{ORDER}_{MESH_TYPE}_tensorial_solver_model{MODEL_TYPE}.cpp
                set(GEN_FILENAME "SEMQ${ORDER}_${MESH_TYPE}_tensorial_solver_${PHYSIC_TYPE}_model${MODEL_TYPE}.cpp")
                set(GEN_FILEPATH "${GENERATED_DIR}/${GEN_FILENAME}")

                # Determine model header and class based on mesh type
                if(MESH_TYPE STREQUAL "struct")
                    set(MODEL_HEADER "model_struct.h")
                    set(MODEL_CLASS "ModelStruct")
                else()
                    set(MODEL_HEADER "model_unstruct.h")
                    set(MODEL_CLASS "ModelUnstruct")
                endif()

                # Determine if this is "onelements" (false) or "onnodes" (true) for last template param
                if(MODEL_TYPE STREQUAL "onelements")
                    set(MODEL_ON_NODES "false")
                else()
                    set(MODEL_ON_NODES "true")
                endif()



                # Determine if ModelStruct needs ORDER template parameter
                if(MESH_TYPE STREQUAL "struct")
                    set(MODEL_TEMPLATE_PARAMS "float, int, ORDER")
                else()
                    set(MODEL_TEMPLATE_PARAMS "float, int")
                endif()

                # Configure the file from template
                configure_file(
                    "${TEMPLATE_FILE}"
                    "${GEN_FILEPATH}"
                    @ONLY
                )

                # Add to list of sources
                list(APPEND TENSORIAL_GENERATED_SOURCES "${GEN_FILEPATH}")

                message(STATUS "Configured: ${GEN_FILENAME}")
            endforeach()
        endforeach()
    endforeach()
endforeach()

set(TENSORIAL_GENERATED_SOURCES ${TENSORIAL_GENERATED_SOURCES} PARENT_SCOPE)
