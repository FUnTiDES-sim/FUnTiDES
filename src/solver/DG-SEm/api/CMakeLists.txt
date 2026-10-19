# Public class/data headers of the DG-SEm solver variant (see DG/api).
target_include_directories(fun_solver
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install_funtides_headers(solver ${CMAKE_CURRENT_SOURCE_DIR}/include)
