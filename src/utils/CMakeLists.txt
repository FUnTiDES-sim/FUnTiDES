add_library(proxy_utils INTERFACE)

    target_include_directories(proxy_utils 
        INTERFACE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/utils/include>
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/discretization/fe/SEMKernels/src>

    )
