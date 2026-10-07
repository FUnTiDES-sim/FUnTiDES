add_library(fun_gradient)
add_library(funtides::gradient ALIAS fun_gradient)
set_target_properties(fun_gradient PROPERTIES EXPORT_NAME gradient)

target_include_directories(fun_gradient PUBLIC
    $<INSTALL_INTERFACE:include/funtides/gradient>
)

add_subdirectory(acoustic)
add_subdirectory(common)
add_subdirectory(elastic)

set(INV_LIBS fun_discretization fun_core fun_model fun_gradient_api)

target_link_kokkos_if_enabled(fun_gradient)

target_link_libraries(fun_gradient PUBLIC ${INV_LIBS})
