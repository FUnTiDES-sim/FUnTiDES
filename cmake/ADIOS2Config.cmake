# ADIOS2.cmake - Configuration for ADIOS2 with minimal dependencies
# All external dependencies disabled for submodule-only build

if(NOT TARGET adios2::adios2)
    message(STATUS "Configuring ADIOS2 with minimal dependencies")

    # Set ADIOS2 build options - disable all external dependencies
    set(ADIOS2_USE_MPI OFF CACHE BOOL "Disable MPI support")
    set(ADIOS2_USE_HDF5 OFF CACHE BOOL "Disable HDF5 engine")
    set(ADIOS2_USE_Python OFF CACHE BOOL "Disable Python bindings")
    set(ADIOS2_USE_Fortran OFF CACHE BOOL "Disable Fortran bindings")
    set(ADIOS2_USE_SST OFF CACHE BOOL "Disable SST engine")
    set(ADIOS2_USE_DataMan OFF CACHE BOOL "Disable DataMan engine")
    set(ADIOS2_USE_DataSpaces OFF CACHE BOOL "Disable DataSpaces engine")
    set(ADIOS2_USE_ZeroMQ OFF CACHE BOOL "Disable ZeroMQ")
    set(ADIOS2_USE_Profiling OFF CACHE BOOL "Disable profiling")
    set(ADIOS2_USE_SysVShMem OFF CACHE BOOL "Disable SysV shared memory")
    set(ADIOS2_USE_UCX OFF CACHE BOOL "Disable UCX")
    set(ADIOS2_USE_IME OFF CACHE BOOL "Disable IME")
    set(ADIOS2_USE_DAOS OFF CACHE BOOL "Disable DAOS")
    set(ADIOS2_USE_GPU_Support OFF CACHE BOOL "Disable GPU support")
    set(ADIOS2_USE_CUDA OFF CACHE BOOL "Disable CUDA")
    set(ADIOS2_USE_Kokkos OFF CACHE BOOL "Disable Kokkos")
    set(ADIOS2_USE_OpenMP OFF CACHE BOOL "Enable OpenMP for compression")
    set(ADIOS2_USE_Blosc2 OFF CACHE BOOL "Enable Blosc compression")
    set(ADIOS2_USE_BZip2 OFF CACHE BOOL "Disable BZip2 compression")
    set(ADIOS2_USE_ZFP OFF CACHE BOOL "Disable ZFP compression")
    set(ADIOS2_USE_SZ OFF CACHE BOOL "Disable SZ compression")
    set(ADIOS2_USE_MGARD OFF CACHE BOOL "Disable MGARD compression")
    set(ADIOS2_USE_PNG OFF CACHE BOOL "Disable PNG support")
    set(ADIOS2_USE_FFTW OFF CACHE BOOL "Disable FFTW")
    set(ADIOS2_USE_Catalyst OFF CACHE BOOL "Disable Catalyst")
    set(ADIOS2_USE_VTK OFF CACHE BOOL "Disable VTK")

    # Build configuration
    set(BUILD_TESTING OFF CACHE BOOL "Disable testing")
    set(ADIOS2_BUILD_TESTING OFF CACHE BOOL "Disable ADIOS2 testing")
    set(ADIOS2_BUILD_EXAMPLES OFF CACHE BOOL "Disable ADIOS2 examples")
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries")

    # Force disable all test building
    set(ADIOS2_BUILD_TESTING_SERIAL OFF CACHE BOOL "Disable serial tests")
    set(ADIOS2_BUILD_TESTING_MPI OFF CACHE BOOL "Disable MPI tests")

    # Prevent CTest from being included
    set(BUILD_TESTING OFF)
    set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)

    # Add ADIOS2 subdirectory
    if(EXISTS ${CMAKE_SOURCE_DIR}/external/ADIOS2/CMakeLists.txt)
        add_subdirectory(${CMAKE_SOURCE_DIR}/external/ADIOS2
                        ${CMAKE_BINARY_DIR}/external/ADIOS2)

    # Create alias for consistent naming
    if(TARGET adios2)
        add_library(adios2::adios2 ALIAS adios2)
    endif()

        message(STATUS "ADIOS2 configured successfully with minimal dependencies")
    else()
        message(WARNING "ADIOS2 submodule not found at external/ADIOS2")
    endif()

endif()
