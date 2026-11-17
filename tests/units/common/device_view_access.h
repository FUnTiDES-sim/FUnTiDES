#ifndef FUNTIDES_TESTS_UNITS_COMMON_DEVICE_VIEW_ACCESS_H_
#define FUNTIDES_TESTS_UNITS_COMMON_DEVICE_VIEW_ACCESS_H_

#include <Kokkos_Core.hpp>
#include <initializer_list>
#include <stdexcept>
#include <string>

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

/**
 * @brief Host-readable view of @p view, to read many entries at once.
 *
 * For a view already in host memory this is @p view itself, not a copy: it follows later changes
 * to @p view. Use copyToHost to keep values across such changes.
 */
template <typename ViewType>
auto toHost(const ViewType& view) {
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
}

/// Host copy of a view that never aliases it, whatever the memory space.
template <typename ViewType>
auto copyToHost(const ViewType& view) {
  auto host = Kokkos::create_mirror(Kokkos::HostSpace{}, view);
  Kokkos::deep_copy(host, view);
  return host;
}

/**
 * @brief Fills a view from the host.
 * @param view View to fill.
 * @param f Called once with a host mirror of @p view, to write into. The mirror holds the current
 *          values of @p view, so entries that @p f leaves untouched keep their value.
 */
template <typename ViewType, typename F>
void fillFromHost(const ViewType& view, F&& f) {
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
  f(host);
  Kokkos::deep_copy(view, host);
}

/// Fills a rank-1 view with the given values, in order; their count must match the view extent.
template <typename ViewType>
void fillFromList(const ViewType& view, std::initializer_list<typename ViewType::non_const_value_type> values) {
  if (values.size() != view.extent(0)) {
    throw std::invalid_argument("fillFromList: " + std::to_string(values.size()) + " values for a view of extent " +
                                std::to_string(view.extent(0)));
  }
  fillFromHost(view, [&](const auto& host) {
    int i = 0;
    for (auto v : values) host(i++) = v;
  });
}

}  // namespace test_utils

#endif  // FUNTIDES_TESTS_UNITS_COMMON_DEVICE_VIEW_ACCESS_H_
