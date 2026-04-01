#include "proxy_macros.h"

#pragma once

// Minimal data_type definition for standalone builds
#ifndef USE_DOUBLE
using real_t = float;
#else
using real_t = double;
#endif