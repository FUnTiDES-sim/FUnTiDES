#ifndef SEM_MACROS_HPP_
#define SEM_MACROS_HPP_

#include "common_config.h"

// define Macros test SEM proxy 2D case
#if defined(SEM2D)
#define DIMENSION 2
#define ROW 36
#define COL 4
#define ZEROED2D 0
// test SEM proxy 3D case
#else
#define DIMENSION 3
#define ROW 64
#define COL 6
#define ZEROED2D 1
#endif

#if defined(USE_KOKKOS) && !defined(SEM_MESHCOLOR)
#define ATOMICADD(ADD1, ADD2) Kokkos::atomic_add(&ADD1, ADD2)
#else
#define ATOMICADD(ADD1, ADD2) ADD1 += ADD2
#endif

// define fence
#if defined(USE_KOKKOS)
#define FENCE Kokkos::fence();
#else
#define FENCE
#endif

#endif  // SEM_MACROS_HPP_
