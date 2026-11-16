# SEM method tree: the physics-generic skeleton in common/ and one
# subdirectory per physics, holding that physics' SEM types and impl
# fragment beside its explicit instantiations. The library target is
# declared in impl/CMakeLists.txt.
add_subdirectory(common)
add_subdirectory(acoustic)
add_subdirectory(elastic)
add_subdirectory(acoustoelastic)
