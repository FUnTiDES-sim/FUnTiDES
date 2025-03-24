#ifndef ARGSPARSE_HPP
#define ARGSPARSE_HPP

#include <string>

char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);

#endif
