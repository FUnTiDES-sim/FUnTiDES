# Boilerplate shared by every physics lives in add_solver_instantiations().
# The p-adaptive solver couples two sub-domains of different polynomial order
# (ORDER_MIN, ORDER_MAX): ORDER_PAIRS dispatches to the generator that emits
# one file per ordered pair instead of one per order.
add_solver_instantiations(
    PHYSIC acoustic
    MAX_ORDER TARGET_SOLVER_ACOUSTIC_ORDER
    FILE_PREFIX DGPADAPTIVEQ
    ORDER_PAIRS
)
