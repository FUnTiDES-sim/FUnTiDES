#ifndef FUNTIDES_UTILS_INCLUDE_DATA_TYPE_KOKKOS_H_
#define FUNTIDES_UTILS_INCLUDE_DATA_TYPE_KOKKOS_H_

#include "Kokkos_Core_fwd.hpp"
#include <Kokkos_Core.hpp>

using MemSpace = Kokkos::SharedSpace;
using Layout = Kokkos::DefaultExecutionSpace::array_layout;

typedef Kokkos::View<int *, Layout, MemSpace> vectorInt;
typedef Kokkos::View<float *, Layout, MemSpace> vectorReal;
typedef Kokkos::View<double *, Layout, MemSpace> vectorDouble;
typedef Kokkos::View<int **, Layout, MemSpace> arrayInt;
typedef Kokkos::View<float **, Layout, MemSpace> arrayReal;
typedef Kokkos::View<double **, Layout, MemSpace> arrayDouble;
typedef Kokkos::View<int ***, Layout, MemSpace> array3DInt;
typedef Kokkos::View<float ***, Layout, MemSpace> array3DReal;
typedef Kokkos::View<double ***, Layout, MemSpace> array3DDouble;
#endif  // FUNTIDES_UTILS_INCLUDE_DATA_TYPE_KOKKOS_H_
