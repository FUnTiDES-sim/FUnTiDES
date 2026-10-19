# SEM solver tree: public API boundary, per-physics implementations, pybind11
# wrapping. All variants compile into the fun_solver library target, which is
# declared in SEM/impl.
add_subdirectory(api)
add_subdirectory(impl)
add_subdirectory(pywrap)
