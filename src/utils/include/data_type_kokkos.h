#ifndef FUNTIDES_UTILS_INCLUDE_DATA_TYPE_KOKKOS_H_
#define FUNTIDES_UTILS_INCLUDE_DATA_TYPE_KOKKOS_H_

#include <Kokkos_Core.hpp>

#include "Kokkos_Core_fwd.hpp"

using Layout = Kokkos::DefaultExecutionSpace::array_layout;
using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;

typedef Kokkos::View<int *, Layout, DeviceSpace> vectorInt;
typedef Kokkos::View<float *, Layout, DeviceSpace> vectorReal;
typedef Kokkos::View<double *, Layout, DeviceSpace> vectorDouble;
typedef Kokkos::View<int **, Layout, DeviceSpace> arrayInt;
typedef Kokkos::View<float **, Layout, DeviceSpace> arrayReal;
typedef Kokkos::View<double **, Layout, DeviceSpace> arrayDouble;
typedef Kokkos::View<int ***, Layout, DeviceSpace> array3DInt;
typedef Kokkos::View<float ***, Layout, DeviceSpace> array3DReal;
typedef Kokkos::View<double ***, Layout, DeviceSpace> array3DDouble;

#endif  // FUNTIDES_UTILS_INCLUDE_DATA_TYPE_KOKKOS_H_
