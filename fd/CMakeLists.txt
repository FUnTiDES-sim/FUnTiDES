################################
# Setup for FD_EXECUTABLE
################################
################################
###############################
###############################

add_executable(fd_executable FDTDmain.cpp)
target_include_directories(fd_executable PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/utilsFDTD)

set(FD_LIBS utils)

if(USE_KOKKOS)
  list(APPEND FD_LIBS Kokkos::kokkos)
endif()

target_link_libraries(fd_executable PRIVATE ${FD_LIBS})
