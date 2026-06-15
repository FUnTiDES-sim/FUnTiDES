add_library(fun_model "")
add_library(funtides::model ALIAS fun_model)

target_link_libraries(fun_model PUBLIC fun_utils Kokkos::kokkos)

add_subdirectory(mesh)
add_subdirectory(grid)
