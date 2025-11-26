add_library(proxy_model_common INTERFACE)

target_include_directories(proxy_model_common
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS proxy_model_common
  EXPORT proxy_model_common_targets
  INCLUDES DESTINATION include)

install(DIRECTORY include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")
