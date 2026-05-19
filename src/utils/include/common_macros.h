#ifndef FUNTIDES_UTILS_INCLUDE_COMMON_MACROS_H_
#define FUNTIDES_UTILS_INCLUDE_COMMON_MACROS_H_

#include "common_config.h"

#define PROXY_HOST_DEVICE KOKKOS_FORCEINLINE_FUNCTION

// MAINLOOP
#define LaunchMaxThreadsPerBlock 128
#define LaunchMinBlocksPerSM 4

// Backward-compatible loop helpers used by legacy FE kernels.
#define LOOPHEAD(Range, Index)         \
  Kokkos::parallel_for(                                                       \
      Kokkos::RangePolicy<>(0, Range),                                        \
      KOKKOS_CLASS_LAMBDA(const int Index)
#define LOOPEND );

#define MAINLOOPHEAD(Range, Index)                                       \
  Kokkos::parallel_for(                                                       \
      "MainLoop",                                                            \
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,      \
                                               LaunchMinBlocksPerSM>>(         \
          0, Range),                                                          \
      KOKKOS_CLASS_LAMBDA(const int Index)
#define MAINLOOPEND );

#define FIND_MAX_1D(Array, Range, Result)                                                          \
  if (Array.extent(0) == 0) throw std::runtime_error("Error in FIND_MAX_1D: Array has zero size"); \
  Kokkos::parallel_reduce(                                                                         \
      "FindMax1D", Range,                                                                          \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_max) {                              \
        if (Array[i] > local_max) local_max = Array[i];                                            \
      },                                                                                           \
      Kokkos::Max<decltype(Result)>(Result));

#define FIND_MIN(Array, Range, Result)                                \
  Kokkos::parallel_reduce(                                            \
      Range,                                                          \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_min) { \
        if (Array[i] < local_min) local_min = Array[i];               \
      },                                                              \
      Kokkos::Min<decltype(Result)>(Result));

#define SUM(Array, Range, Result)                                                                      \
  Kokkos::parallel_reduce(                                                                             \
      Range, KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_sum) { local_sum += Array[i]; }, \
      Kokkos::Sum<decltype(Result)>(Result));

#define KOKKOSNAME "v",

#endif  // FUNTIDES_UTILS_INCLUDE_COMMON_MACROS_H_
