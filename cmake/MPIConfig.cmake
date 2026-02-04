if(USE_MPI)
find_package(MPI REQUIRED)
add_compile_definitions(USE_MPI)
endif()
