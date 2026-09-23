#ifndef FUNTIDES_CORE_INCLUDE_COMMON_MACROS_H_
#define FUNTIDES_CORE_INCLUDE_COMMON_MACROS_H_

#include "common_config.h"

/// Qualifier for functions callable from host and device code.
#define PROXY_HOST_DEVICE KOKKOS_FORCEINLINE_FUNCTION

/// Kokkos launch bound: maximum threads per block for MAINLOOPHEAD kernels.
#define LaunchMaxThreadsPerBlock 128
/// Kokkos launch bound: minimum resident blocks per SM for MAINLOOPHEAD kernels.
#define LaunchMinBlocksPerSM 4

/**
 * @brief Opens a Kokkos::parallel_for over [0, Range) on the default execution space.
 *
 * The macro leaves the call open: the loop body (a block) follows and must be closed
 * with LOOPEND. The body is a KOKKOS_CLASS_LAMBDA, so it captures `this`.
 * @param Range Number of iterations.
 * @param Index Name of the `const int` loop index available in the body.
 */
#define LOOPHEAD(Range, Index)         \
  Kokkos::parallel_for(                                                       \
      Kokkos::RangePolicy<>(0, Range),                                        \
      KOKKOS_CLASS_LAMBDA(const int Index)
/// Closes a loop opened with LOOPHEAD.
#define LOOPEND );

/**
 * @brief Opens a named Kokkos::parallel_for ("MainLoop") over [0, Range) with launch bounds
 * LaunchMaxThreadsPerBlock and LaunchMinBlocksPerSM.
 *
 * The body follows and must be closed with MAINLOOPEND. The body is a KOKKOS_CLASS_LAMBDA,
 * so it captures `this`.
 * @param Range Number of iterations.
 * @param Index Name of the `const int` loop index available in the body.
 */
#define MAINLOOPHEAD(Range, Index)                                       \
  Kokkos::parallel_for(                                                       \
      "MainLoop",                                                            \
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,      \
                                               LaunchMinBlocksPerSM>>(         \
          0, Range),                                                          \
      KOKKOS_CLASS_LAMBDA(const int Index)
/// Closes a loop opened with MAINLOOPHEAD.
#define MAINLOOPEND );

/**
 * @brief Computes the maximum of a 1D array with a Kokkos::parallel_reduce.
 *
 * Expands to several statements and uses KOKKOS_CLASS_LAMBDA (captures `this`).
 * @param Array Indexable array-like object with extent(0).
 * @param Range Number of elements to scan, from index 0.
 * @param Result Scalar lvalue receiving the maximum; its type sets the reduction type.
 * @throws std::runtime_error If Array has zero extent.
 */
#define FIND_MAX_1D(Array, Range, Result)                                                          \
  if (Array.extent(0) == 0) throw std::runtime_error("Error in FIND_MAX_1D: Array has zero size"); \
  Kokkos::parallel_reduce(                                                                         \
      "FindMax1D", Range,                                                                          \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_max) {                              \
        if (Array[i] > local_max) local_max = Array[i];                                            \
      },                                                                                           \
      Kokkos::Max<decltype(Result)>(Result));

/**
 * @brief Computes the minimum of a 1D array with a Kokkos::parallel_reduce.
 *
 * Uses KOKKOS_CLASS_LAMBDA (captures `this`).
 * @param Array Indexable array-like object.
 * @param Range Number of elements to scan, from index 0.
 * @param Result Scalar lvalue receiving the minimum; its type sets the reduction type.
 */
#define FIND_MIN(Array, Range, Result)                                \
  Kokkos::parallel_reduce(                                            \
      Range,                                                          \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_min) { \
        if (Array[i] < local_min) local_min = Array[i];               \
      },                                                              \
      Kokkos::Min<decltype(Result)>(Result));

/**
 * @brief Computes the sum of a 1D array with a Kokkos::parallel_reduce.
 *
 * Uses KOKKOS_CLASS_LAMBDA (captures `this`).
 * @param Array Indexable array-like object.
 * @param Range Number of elements to sum, from index 0.
 * @param Result Scalar lvalue receiving the sum; its type sets the reduction type.
 */
#define SUM(Array, Range, Result)                                                                      \
  Kokkos::parallel_reduce(                                                                             \
      Range, KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_sum) { local_sum += Array[i]; }, \
      Kokkos::Sum<decltype(Result)>(Result));

/// @todo VERIFY: what is KOKKOSNAME used for? It expands to the string literal "v" followed by a comma, apparently as the label argument of a Kokkos::View constructor.
#define KOKKOSNAME "v",

#endif  // FUNTIDES_CORE_INCLUDE_COMMON_MACROS_H_
