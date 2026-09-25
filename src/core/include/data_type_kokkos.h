#ifndef FUNTIDES_CORE_INCLUDE_DATA_TYPE_KOKKOS_H_
#define FUNTIDES_CORE_INCLUDE_DATA_TYPE_KOKKOS_H_

#include <Kokkos_Core.hpp>

#include "Kokkos_Core_fwd.hpp"

/// Array layout of the default Kokkos execution space.
using Layout = Kokkos::DefaultExecutionSpace::array_layout;
/// Memory space of the default Kokkos execution space (device memory on GPU builds).
using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

/// @name Kokkos views in DeviceSpace with the default Layout
/// Rank 1, 2 and 3 views of int, float and double.
/// @{
typedef Kokkos::View<int *, Layout, DeviceSpace> vectorInt;
typedef Kokkos::View<float *, Layout, DeviceSpace> vectorReal;
typedef Kokkos::View<double *, Layout, DeviceSpace> vectorDouble;
typedef Kokkos::View<int **, Layout, DeviceSpace> arrayInt;
typedef Kokkos::View<float **, Layout, DeviceSpace> arrayReal;
typedef Kokkos::View<double **, Layout, DeviceSpace> arrayDouble;
typedef Kokkos::View<int ***, Layout, DeviceSpace> array3DInt;
typedef Kokkos::View<float ***, Layout, DeviceSpace> array3DReal;
typedef Kokkos::View<double ***, Layout, DeviceSpace> array3DDouble;
/// @}

#endif  // FUNTIDES_CORE_INCLUDE_DATA_TYPE_KOKKOS_H_
