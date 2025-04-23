#ifndef ARGSPARSE_HPP
#define ARGSPARSE_HPP

#include <string>
#include <algorithm> // Pour std::find

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
inline char *getCmdOption(char **begin, char **end, const std::string &option) {
  char **itr = std::find(begin, end, option.c_str());
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
inline bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option.c_str()) != end;
}

#endif // ARGSPARSE_HPP
