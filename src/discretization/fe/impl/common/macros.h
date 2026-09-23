/**
 * @file
 * @brief Function qualifiers for the math helpers of mathUtilites.h.
 *
 * SEMKERNELS_HOST_DEVICE makes a function callable from host and device code
 * when compiled by nvcc or hipcc, and expands to nothing otherwise.
 * SEMKERNELS_INLINE expands to `inline`.
 */
#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define SEMKERNELS_HOST_DEVICE __host__ __device__
#else
#define SEMKERNELS_HOST_DEVICE
#endif

#define SEMKERNELS_INLINE inline
