#pragma once


#if defined(__CUDACC__) || defined(__HIPCC__)
#define SEMKERNELS_HOST_DEVICE __host__ __device__
#else
#define SEMKERNELS_HOST_DEVICE
#endif

#define SEMKERNELS_INLINE inline
