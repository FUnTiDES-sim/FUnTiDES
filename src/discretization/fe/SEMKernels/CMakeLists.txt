cmake_minimum_required(VERSION 3.22)

set(SEM_KERNELS_VERSION_MAJOR 0)
set(SEM_KERNELS_VERSION_MINOR 1)
set(SEM_KERNELS_VERSION_PATCHLEVEL 0)


# Detect if submodule or standalone
get_directory_property(parent_dir PARENT_DIRECTORY)
if(parent_dir)
    set(SEMKERNELS_IS_SUBMODULE TRUE)
    option(SEMKERNELS_BUILD_TESTS "Build SEMKernels tests" OFF)
else()
    set(SEMKERNELS_IS_SUBMODULE FALSE)
    add_library(discretization INTERFACE)
    project(SEMKernels VERSION ${SEM_KERNELS_VERSION_MAJOR}.${SEM_KERNELS_VERSION_MINOR}.${SEM_KERNELS_VERSION_PATCHLEVEL} LANGUAGES CXX)
    option(SEMKERNELS_BUILD_TESTS "Build SEMKernels tests" ON)
endif()

# Configure includes based on context
if(SEMKERNELS_IS_SUBMODULE)
    message(STATUS "SEMKernels: Using FUnTiDES headers (submodule mode)")
else()
    message(STATUS "SEMKernels: Using standalone external headers")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)
endif()

# Always build the library
add_subdirectory(src)

# Build tests only if requested
if(SEMKERNELS_BUILD_TESTS)
    # Setup GTest
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        message(STATUS "GTest not found, fetching from GitHub...")
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()
    
    include(GoogleTest)
    
    if(NOT COMMAND add_gtest)
        function(add_gtest TEST_NAME TEST_FILE EXTRA_LINKS)
            add_executable(${TEST_NAME} ${TEST_FILE})
            target_link_libraries(${TEST_NAME} PRIVATE GTest::gtest GTest::gtest_main ${EXTRA_LINKS})
            gtest_discover_tests(${TEST_NAME})
        endfunction()
    endif()
    
    enable_testing()
    add_subdirectory(tests)
endif()

if(SEMKERNELS_ENABLE_DOCS)
    # add_subdirectory(docs)
endif()