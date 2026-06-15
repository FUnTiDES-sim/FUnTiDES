# ============================================================================
# Include the generator function
# ============================================================================
include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

# ============================================================================
# Configuration: Define all combinations to generate
# ============================================================================

# Define polynomial orders
set(ORDERS "")
foreach(order RANGE 1 ${TARGET_SOLVER_ACOUSTIC_ORDER})
    list(APPEND ORDERS ${order})
endforeach()

# Define mesh types
set(MESH_TYPES "struct" "unstruct")

# Define model types
set(MODEL_TYPES "onelements" "onnodes")

# ============================================================================
# Generate solver implementation files
# ============================================================================

generate_solver_implementations(
    SOLVER_NAME "makutu"
    ORDERS ${ORDERS}
    MESH_TYPES ${MESH_TYPES}
    MODEL_TYPES ${MODEL_TYPES}
    PHYSIC_TYPES "acoustic"
    OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    OUTPUT_VAR MAKUTU_GENERATED_SOURCES
)

target_sources(fun_solver PRIVATE ${MAKUTU_GENERATED_SOURCES})
