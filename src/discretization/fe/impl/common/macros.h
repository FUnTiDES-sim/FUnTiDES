/**
 * @file
 * @brief Host/device and inline qualifier macros for the math helpers of mathUtilites.h.
 */
#pragma once

/// Makes a function callable from host and device code under nvcc or hipcc; empty otherwise.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define SEMKERNELS_HOST_DEVICE __host__ __device__
#else
#define SEMKERNELS_HOST_DEVICE
#endif

/// Expands to `inline`.
#define SEMKERNELS_INLINE inline
