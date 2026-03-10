#ifndef FUNTIDES_SOLVER_FD_INCLUDE_FD_SOURCE_RECEIVERS_H_
#define FUNTIDES_SOLVER_FD_INCLUDE_FD_SOURCE_RECEIVERS_H_
#include <data_type.h>

#include "fd_macros.h"
#include "fd_options.h"

struct FdtdSourceReceivers
{
  // source location
  int xsrc{-1}, ysrc{-1}, zsrc{-1};
};
#endif  // FUNTIDES_SOLVER_FD_INCLUDE_FD_SOURCE_RECEIVERS_H_
