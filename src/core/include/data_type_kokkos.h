#ifndef FUNTIDES_CORE_INCLUDE_DATA_TYPE_KOKKOS_H_
#define FUNTIDES_CORE_INCLUDE_DATA_TYPE_KOKKOS_H_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <clocale>

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

/// @name Kokkos DualViews: one copy in DeviceSpace, one in host memory
///
/// view_device() has exactly the type of the matching view above (dualVectorReal::t_dev is
/// vectorReal), so it can be passed wherever a device view is expected.
///
/// Protocol: write one side, declare it with modify_host() / modify_device(), then call
/// sync_device() / sync_host() before the other side is read. A sync only copies when the other
/// side was declared modified, so calling it when nothing changed costs nothing. Without a device
/// backend both sides share the same memory and these calls compile to nothing.
/// @{
typedef Kokkos::DualView<int *, Layout, DeviceSpace> dualVectorInt;
typedef Kokkos::DualView<float *, Layout, DeviceSpace> dualVectorReal;
typedef Kokkos::DualView<double *, Layout, DeviceSpace> dualVectorDouble;
typedef Kokkos::DualView<int **, Layout, DeviceSpace> dualArrayInt;
typedef Kokkos::DualView<float **, Layout, DeviceSpace> dualArrayReal;
typedef Kokkos::DualView<double **, Layout, DeviceSpace> dualArrayDouble;
typedef Kokkos::DualView<int ***, Layout, DeviceSpace> dualArray3DInt;
typedef Kokkos::DualView<float ***, Layout, DeviceSpace> dualArray3DReal;
typedef Kokkos::DualView<double ***, Layout, DeviceSpace> dualArray3DDouble;
/// @}

/**
 * @brief Reads one entry of a DualView from the side the caller runs on.
 *
 * In device code it reads view_device(), in host code view_host(), so the caller never touches
 * memory it cannot access. The side that is read must be up to date: after a modification on
 * the other side, call sync_host() or sync_device() first.
 *
 * @tparam DualViewType Kokkos::DualView type.
 * @tparam Indices Integral index types, one per rank.
 * @param[in] dv DualView to read.
 * @param[in] indices Indices of the entry.
 * @return Value of the entry.
 */
template <class DualViewType, class... Indices>
KOKKOS_FORCEINLINE_FUNCTION typename DualViewType::t_dev::non_const_value_type dualAt(const DualViewType &dv,
                                                                                      const Indices... indices) {
  KOKKOS_IF_ON_DEVICE((return dv.view_device()(indices...);))
  KOKKOS_IF_ON_HOST((return dv.view_host()(indices...);))
}

/**
 * @brief Writes one entry of a DualView on the side the caller runs on.
 *
 * In device code it writes view_device(), in host code view_host(). It does not mark the side as
 * modified: the caller does it once, from the host, after the writes (modify_device() after a
 * kernel, modify_host() after a host loop).
 *
 * @tparam DualViewType Kokkos::DualView type.
 * @tparam Indices Integral index types, one per rank.
 * @param[in] dv DualView to write (the views it holds are shared, hence the const reference).
 * @param[in] value Value to store.
 * @param[in] indices Indices of the entry.
 */
template <class DualViewType, class... Indices>
KOKKOS_FORCEINLINE_FUNCTION void dualSet(const DualViewType &dv,
                                         const typename DualViewType::t_dev::non_const_value_type value,
                                         const Indices... indices) {
  KOKKOS_IF_ON_DEVICE((dv.view_device()(indices...) = value;))
  KOKKOS_IF_ON_HOST((dv.view_host()(indices...) = value;))
}

#endif  // FUNTIDES_CORE_INCLUDE_DATA_TYPE_KOKKOS_H_
