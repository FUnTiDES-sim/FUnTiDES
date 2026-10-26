#ifndef FUNTIDES_TESTS_UNITS_COMMON_DEVICE_VIEW_ACCESS_H_
#define FUNTIDES_TESTS_UNITS_COMMON_DEVICE_VIEW_ACCESS_H_

#include <Kokkos_Core.hpp>
#include <initializer_list>

/**
 * @file
 * @brief Host access to the entries of device views, for unit tests.
 *
 * Views allocated in a device memory space cannot be dereferenced on the host. These helpers go
 * through Kokkos::deep_copy instead; with a host memory space they cost no copy.
 */

namespace test_utils {

/// Reads one entry of a view, wherever it lives.
template <typename ViewType, typename... Indices>
typename ViewType::non_const_value_type getAt(const ViewType& view, Indices... indices) {
  typename ViewType::non_const_value_type value{};
  Kokkos::deep_copy(value, Kokkos::subview(view, indices...));
  return value;
}

/// Writes one entry of a view, wherever it lives.
template <typename ViewType, typename ValueType, typename... Indices>
void setAt(const ViewType& view, ValueType value, Indices... indices) {
  Kokkos::deep_copy(Kokkos::subview(view, indices...), static_cast<typename ViewType::non_const_value_type>(value));
}

/// Host copy of a view, to read many entries at once.
template <typename ViewType>
auto toHost(const ViewType& view) {
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
}

/**
 * @brief Fills a view from the host.
 * @param view View to fill.
 * @param f Called once with a host mirror of @p view, to write into.
 */
template <typename ViewType, typename F>
void fillFromHost(const ViewType& view, F&& f) {
  auto host = Kokkos::create_mirror_view(view);
  f(host);
  Kokkos::deep_copy(view, host);
}

/// Fills a rank-1 view with the given values, in order.
template <typename ViewType>
void fillFromList(const ViewType& view, std::initializer_list<typename ViewType::non_const_value_type> values) {
  fillFromHost(view, [&](const auto& host) {
    int i = 0;
    for (auto v : values) host(i++) = v;
  });
}

}  // namespace test_utils

#endif  // FUNTIDES_TESTS_UNITS_COMMON_DEVICE_VIEW_ACCESS_H_
