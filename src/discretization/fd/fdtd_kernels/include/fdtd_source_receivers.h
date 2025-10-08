#ifndef FDTD_SOURCE_RECEIVERS_H_
#define FDTD_SOURCE_RECEIVERS_H_

#include "fdtd_macros.h"
#include "fdtd_options.h"
#include <data_type.h>

struct FdtdSourceReceivers
{
  // source location
  int xsrc{-1},ysrc{-1}, zsrc{-1};
};
#endif //FDTD_SOURCE_RECEIVERS_H_
