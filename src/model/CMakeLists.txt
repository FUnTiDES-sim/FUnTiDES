add_library(fun_model "")
add_library(funtides::model ALIAS fun_model)
set_target_properties(fun_model PROPERTIES EXPORT_NAME model)

target_link_libraries(fun_model PUBLIC fun_core Kokkos::kokkos)

target_include_directories(fun_model PUBLIC
    $<INSTALL_INTERFACE:include/funtides/model>
)

add_subdirectory(mesh)
