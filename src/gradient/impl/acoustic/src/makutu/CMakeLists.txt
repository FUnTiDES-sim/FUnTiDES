# ============================================================================
# Include the generator function
# ============================================================================
include(${CMAKE_SOURCE_DIR}/cmake/GenerateDifferentiatorImplementations.cmake)

# ============================================================================
# Configuration: Define all combinations to generate
# ============================================================================

set(ORDERS "")
foreach(order RANGE 1 ${TARGET_DIFFERENTIATOR_ACOUSTIC_ORDER})
    list(APPEND ORDERS ${order})
endforeach()

# Define mesh types
set(MESH_TYPES "struct" "unstruct")

# Define model types
set(MODEL_TYPES "onelements" "onnodes")

# Define physics types
set(PHYSIC_TYPES "acoustic")

# ============================================================================
# Generate differentiator implementation files
# ============================================================================

generate_differentiator_implementations(
    ORDERS ${ORDERS}
    MESH_TYPES ${MESH_TYPES}
    MODEL_TYPES ${MODEL_TYPES}
    PHYSIC_TYPES ${PHYSIC_TYPES}
    OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    OUTPUT_VAR DIFFERENTIATOR_MAKUTU_GENERATED_SOURCES
)

target_sources(funtides_gradient PRIVATE ${DIFFERENTIATOR_MAKUTU_GENERATED_SOURCES})
