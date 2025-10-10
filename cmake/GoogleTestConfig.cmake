#-------------------------------------------------------------------
# Google test configuration
#-------------------------------------------------------------------

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/googletest/CMakeLists.txt")
    message(FATAL_ERROR "googletest submodule not found. Run 'git submodule update --init --recursive'")
else()
    add_subdirectory(${CMAKE_SOURCE_DIR}/external/googletest)
    include(GoogleTest)
    enable_testing()
endif()