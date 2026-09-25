#ifndef FUNTIDES_CORE_INCLUDE_SEM_MACROS_H_
#define FUNTIDES_CORE_INCLUDE_SEM_MACROS_H_
#include "common_config.h"

/** @brief Number of spatial dimensions. */
#define DIMENSION 3
/** @brief @todo VERIFY: what does ROW (64) size, and in which arrays is it used? */
#define ROW 64
/** @brief @todo VERIFY: what does COL (6) size, and in which arrays is it used? */
#define COL 6
/** @brief @todo VERIFY: meaning of ZEROED2D (1); is it a flag and who reads it? */
#define ZEROED2D 1

/**
 * @brief Atomic addition on a device or host value.
 * @param ADD1 Lvalue that receives the sum; its address is taken.
 * @param ADD2 Value to add.
 */
#define ATOMICADD(ADD1, ADD2) Kokkos::atomic_add(&ADD1, ADD2)

/** @brief Global Kokkos fence: waits for all pending device work. */
#define FENCE Kokkos::fence();
#endif  // FUNTIDES_CORE_INCLUDE_SEM_MACROS_H_
