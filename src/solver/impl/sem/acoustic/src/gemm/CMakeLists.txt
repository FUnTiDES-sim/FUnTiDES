# Explicit instantiations of the SEM acoustic solver (gemm back-end), for the
# gradient/GEMM proxy path.
add_solver_instantiations(PHYSIC acoustic MAX_ORDER TARGET_SOLVER_ACOUSTIC_ORDER SOLVER_NAME gemm TARGET proxy_solver)
