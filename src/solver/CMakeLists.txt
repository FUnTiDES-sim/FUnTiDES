# The solver library: the pure contract in api/, every implementation under
# impl/ (construction entry point and shared per-physics types in
# impl/common, one subtree per discretization method), and the pyfuntides
# bindings in pywrap/.
add_subdirectory(api)
add_subdirectory(impl)
add_subdirectory(pywrap)
