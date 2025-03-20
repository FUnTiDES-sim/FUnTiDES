//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "SEMproxy.hpp"
#ifdef USE_CALIPER
#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#endif
 
/**
 * @brief Retrieves the value associated with a command-line option.
 * 
 * This function searches for a specified option in the command-line arguments 
 * and returns the corresponding value if found.
 * 
 * @param begin Pointer to the beginning of the argument list.
 * @param end Pointer to the end of the argument list.
 * @param option The option to search for.
 * @return A pointer to the value associated with the option if found, 
 *         otherwise returns nullptr.
 */
char *getCmdOption(char **begin, char **end, const std::string &option) {
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return nullptr;
}

/**
 * @brief Checks if a specific command-line option exists.
 * 
 * This function determines whether a given option is present in the 
 * command-line arguments.
 * 
 * @param begin Pointer to the beginning of the argument list.
 * @param end Pointer to the end of the argument list.
 * @param option The option to check for.
 * @return True if the option exists, otherwise false.
 */
bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[]) {

  time_point<system_clock> startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  cout << "Using Kokkos" << endl;
  Kokkos::initialize(argc, argv);
  {
#endif

#ifdef USE_CALIPER
    cali::ConfigManager mgr;

    char* cali_configuration;
    if (cmdOptionExists(argv, argv + argc, "-P")) {
      cali_configuration = getCmdOption(argv, argc+argv, "-P");
    }
    else {
      cali_configuration = "runtime-report";
    }

    mgr.add(cali_configuration);
    if (mgr.error()) {
      std::cerr << "Config error: " << mgr.error_msg() << std::endl;
    }
    else {
      std::cout << "Starting Caliper with option set at: " << cali_configuration << std::endl;
    }
    mgr.start();
    CALI_CXX_MARK_FUNCTION;
#endif

    cout << "\n+================================= " << endl;
    cout << "| Initializing SEM Application ... " << endl;
    cout << "+================================= \n" << endl;

    SEMproxy semsim(argc, argv);

    semsim.initFiniteElem();

    cout << "\n+================================= " << endl;
    cout << "| Running SEM Application ...      " << endl;
    cout << "+================================= \n" << endl;

    // start timer
    time_point<system_clock> startRunTime = system_clock::now();
    semsim.run();

    cout << "\n+================================= " << endl;
    cout << "| SEM Application Finished.       " << endl;
    cout << "+================================= \n" << endl;

    // print timing information
    cout << "Elapsed Initial Time : "
         << (startRunTime - startInitTime).count() / 1E9 << " seconds." << endl;
    cout << "Elapsed Compute Time : "
         << (system_clock::now() - startRunTime).count() / 1E9 << " seconds."
         << endl;

#ifdef USE_CALIPER
    mgr.flush();
#endif // USE_CALIPER

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExe Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
