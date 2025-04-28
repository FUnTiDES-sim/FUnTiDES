################################
# Setup for FD_EXECUTABLE
################################
################################
###############################
###############################

add_executable(fd_executable FDTDmain.cpp)
target_include_directories(fd_executable PRIVATE ${CMAKE_SOURCE_DIR}/finit_elem/include)
