# ADIOS2.cmake - Minimal Configuration for BP5 File I/O only

if(NOT TARGET adios2::adios2)
    message(STATUS "Configuring ADIOS2 for minimal BP5 File I/O")

    # 1. Disable Parallelism and Language Bindings
    if(USE_MPI)
        message(STATUS "ADIOS2 is configured with MPI")
        set(ADIOS2_USE_MPI ON CACHE BOOL "Enable MPI" FORCE)
    else()
        message(STATUS "ADIOS2 is configured without MPI")
        set(ADIOS2_USE_MPI OFF CACHE BOOL "Disable MPI" FORCE)
    endif()
    set(ADIOS2_USE_Python OFF CACHE BOOL "Disable Python" FORCE)
    set(ADIOS2_USE_Fortran OFF CACHE BOOL "Disable Fortran" FORCE)

    # 2. Disable Network/Staging Engines (Removes EVPath, ENet, Dill, FFS, Atl)
    set(ADIOS2_USE_SST OFF CACHE BOOL "Disable SST Staging" FORCE)
    set(ADIOS2_USE_SSC OFF CACHE BOOL "Disable SSC Staging" FORCE)
    set(ADIOS2_USE_DataMan OFF CACHE BOOL "Disable DataMan" FORCE)
    set(ADIOS2_USE_DataSpaces OFF CACHE BOOL "Disable DataSpaces" FORCE)
    set(ADIOS2_USE_ZeroMQ OFF CACHE BOOL "Disable ZeroMQ" FORCE)
    set(ADIOS2_USE_HDF5 OFF CACHE BOOL "Disable HDF5" FORCE)

    # 3. Disable Legacy/Unused Engines (Keep only BP5)
    set(ADIOS2_USE_BP3 OFF CACHE BOOL "Disable BP3" FORCE)
    set(ADIOS2_USE_BP4 OFF CACHE BOOL "Disable BP4" FORCE)
    # Note: BP5 is usually ON by default, but we leave it implicitly enabled.

    # 4. Disable Compression & External Libs (Speeds up compile)
    set(ADIOS2_USE_Blosc2 OFF CACHE BOOL "Disable Blosc2" FORCE)
    set(ADIOS2_USE_BZip2 OFF CACHE BOOL "Disable BZip2" FORCE)
    set(ADIOS2_USE_ZFP OFF CACHE BOOL "Disable ZFP" FORCE)
    set(ADIOS2_USE_SZ OFF CACHE BOOL "Disable SZ" FORCE)
    set(ADIOS2_USE_PNG OFF CACHE BOOL "Disable PNG" FORCE)
    set(ADIOS2_USE_MGARD OFF CACHE BOOL "Disable MGARD" FORCE)

    # 5. Disable System/Hardware features
    set(ADIOS2_USE_SysVShMem OFF CACHE BOOL "Disable SysV Shared Mem" FORCE)
    set(ADIOS2_USE_UCX OFF CACHE BOOL "Disable UCX" FORCE)
    set(ADIOS2_USE_DAOS OFF CACHE BOOL "Disable DAOS" FORCE)
    set(ADIOS2_USE_IME OFF CACHE BOOL "Disable IME" FORCE)
    set(ADIOS2_USE_CUDA OFF CACHE BOOL "Disable CUDA in ADIOS" FORCE)
    set(ADIOS2_USE_Kokkos OFF CACHE BOOL "Disable Kokkos in ADIOS" FORCE)
    set(ADIOS2_USE_Profiling OFF CACHE BOOL "Disable internal profiling" FORCE)

    # 6. Disable Tools, Tests and Examples
    set(ADIOS2_BUILD_EXAMPLES OFF CACHE BOOL "Disable Examples" FORCE)
    set(ADIOS2_BUILD_TESTING OFF CACHE BOOL "Disable Testing" FORCE)
    set(ADIOS2_INSTALL_GENERATE_CONFIG OFF CACHE BOOL "Disable Config Generation" FORCE)

    # Prevent ADIOS from building its own tools (bpls, adios_reorganize, etc)
    # This saves significant compilation time
    set(ADIOS2_TOOL_BPLS OFF CACHE BOOL "Disable bpls tool" FORCE)
    set(ADIOS2_TOOL_REORGANIZE OFF CACHE BOOL "Disable reorganize tool" FORCE)

    # --------------------------------------------------------
    # Add Subdirectory
    # --------------------------------------------------------
    if(EXISTS ${CMAKE_SOURCE_DIR}/external/ADIOS2/CMakeLists.txt)
        # Prevent CTest from interfering
        set(BUILD_TESTING_SAVED ${BUILD_TESTING})
        set(BUILD_TESTING OFF)

        add_subdirectory(${CMAKE_SOURCE_DIR}/external/ADIOS2
                        ${CMAKE_BINARY_DIR}/external/ADIOS2
                        EXCLUDE_FROM_ALL) # Prevents installing unused artifacts

        set(BUILD_TESTING ${BUILD_TESTING_SAVED})

        if(TARGET adios2_cxx11)
            add_library(adios2::adios2 ALIAS adios2_cxx11)
            message(STATUS "ADIOS2 configured. Target: adios2::adios2")
        endif()
    else()
        message(WARNING "ADIOS2 submodule not found at external/ADIOS2")
    endif()

endif()
