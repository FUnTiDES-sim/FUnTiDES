#ifndef DATATYPE_KOKKOS_HPP_
#define DATATYPE_KOKKOS_HPP_

#include "Kokkos_Core_fwd.hpp"

#ifdef ENABLE_HIP
#define __HIP_PLATFORM_AMD__ 1
#endif

#include <Kokkos_Core.hpp>

#ifdef ENABLE_CUDA
#define MemSpace Kokkos::SharedSpace
using Layout = Kokkos::LayoutLeft;
#else
#define MemSpace Kokkos::HostSpace
using Layout = Kokkos::LayoutRight;
#endif

typedef Kokkos::View<int *, Layout, MemSpace> vectorInt;
typedef Kokkos::View<float *, Layout, MemSpace> vectorReal;
typedef Kokkos::View<double *, Layout, MemSpace> vectorDouble;

typedef Kokkos::View<int **, Layout, MemSpace> arrayInt;
typedef Kokkos::View<float **, Layout, MemSpace> arrayReal;
typedef Kokkos::View<double **, Layout, MemSpace> arrayDouble;

typedef Kokkos::View<int ***, Layout, MemSpace> array3DInt;
typedef Kokkos::View<float ***, Layout, MemSpace> array3DReal;
typedef Kokkos::View<double ***, Layout, MemSpace> array3DDouble;

#endif  // DATATYPE_KOKKOS_HPP_
