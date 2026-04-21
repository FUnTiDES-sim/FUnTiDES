# ============================================================================
# Include the generator function
# ============================================================================
include(${CMAKE_SOURCE_DIR}/cmake/GenerateSolverImplementations.cmake)

# ============================================================================
# Configuration: Define all combinations to generate
# ============================================================================

# Define polynomial orders
set(ORDERS 1 2 3 4 5)

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
    PHYSIC_TYPES "elastic"
    OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    OUTPUT_VAR MAKUTU_GENERATED_SOURCES
)

target_sources(proxy_solver PRIVATE ${MAKUTU_GENERATED_SOURCES})
