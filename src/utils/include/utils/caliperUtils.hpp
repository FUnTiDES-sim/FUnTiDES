#ifndef CALIPER_UTILS_HPP
#define CALIPER_UTILS_HPP

#ifdef USE_CALIPER
#include "argsparse.hpp"
#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#include <string>

inline int launch_caliper_ctx(int argc, char **argv, cali::ConfigManager & mgr) {
  // Defines CALIPER configuration with cmdline arguments
  // Sets -P option to choose caliper outputs
  std::string cali_configuration;
  if (cmdOptionExists(argv, argv + argc, "-P")) {
    cali_configuration = getCmdOption(argv, argc + argv, "-P");
  } else {
    cali_configuration = "runtime-report";
  }

  mgr.add(cali_configuration.c_str());
  if (mgr.error()) {
    std::cerr << "Config error: " << mgr.error_msg() << std::endl;
    return -1;
  } else {
    std::cout << "Starting Caliper with option set at: " << cali_configuration
              << std::endl;
  }
  mgr.start();
  return 0;
}

#endif // USE_CALIPER

#endif // CALIPER_UTILS_HPP
