#ifndef FUNTIDES_UTILS_INCLUDE_SEM_MACROS_H_
#define FUNTIDES_UTILS_INCLUDE_SEM_MACROS_H_
#include "common_config.h"

#define DIMENSION 3
#define ROW 64
#define COL 6
#define ZEROED2D 1

#define ATOMICADD(ADD1, ADD2) Kokkos::atomic_add(&ADD1, ADD2)
/* #define ATOMICADD(ADD1, ADD2) ((ADD1) += (ADD2)) */

#define FENCE Kokkos::fence();
#endif  // FUNTIDES_UTILS_INCLUDE_SEM_MACROS_H_
