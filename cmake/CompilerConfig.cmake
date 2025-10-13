#-------------------------------------------------------------------
# Compilation flags
#-------------------------------------------------------------------

# CC
if(DEFINED CMAKE_C_COMPILER)
  message(STATUS "The CMAKE_C_COMPILER is ${CMAKE_C_COMPILER}")
else()
  message(STATUS "CMAKE_C_COMPILER not set, letting CMake autodetect.")
endif()
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "" FORCE)

# C++
if(DEFINED CMAKE_CXX_COMPILER)
  message(STATUS "The CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
else()
  message(STATUS "CMAKE_CXX_COMPILER not set, letting CMake autodetect.")
endif()

# Add -lstdc++fs if required (not default for some compilers)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-lstdc++fs" HAS_LIBSTDCXXFS)
message(STATUS "Checking for -lstdc++fs: ${HAS_LIBSTDCXXFS}")
if(HAS_LIBSTDCXXFS)
  set(STDCXXFS_LIB stdc++fs) # filesystem support for C++17 (nvhpc, gcc<9)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "" FORCE)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable gnu++ extensions" FORCE)

# CUDA
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "" FORCE)
set(CMAKE_CUDA_ARCHITECTURES native CACHE STRING "Auto-detect CUDA arch" FORCE)