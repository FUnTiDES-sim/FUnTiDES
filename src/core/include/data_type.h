#ifndef FUNTIDES_CORE_INCLUDE_DATA_TYPE_H_
#define FUNTIDES_CORE_INCLUDE_DATA_TYPE_H_
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "common_macros.h"
#include "data_type_kokkos.h"

/// Floating-point type of the solver. Double if USE_DOUBLE is defined, float otherwise.
#ifdef USE_DOUBLE
using real_t = double;
#else
using real_t = float;
#endif

/**
 * @brief Allocates a 1D Kokkos view.
 * @tparam T Kokkos view type of rank 1.
 * @param[in] n1 Number of entries.
 * @return The allocated view.
 */
template <class T>
T allocateVector(int n1) {
  T vect(KOKKOSNAME n1);
  return vect;
}
/**
 * @brief Allocates a 1D Kokkos view or DualView labeled @p name.
 *
 * The label appears in Kokkos error messages (out-of-bounds or inaccessible memory space) and
 * in the Kokkos Tools, so it should identify the array.
 * @tparam T Kokkos view or DualView type of rank 1.
 * @param[in] n1 Number of entries.
 * @param[in] name Label of the allocation.
 * @return The allocated view, zero-initialized.
 */
template <class T>
T allocateVector(int n1, const char *name) {
  T vect(std::string(name), n1);  // View labels must be std::string or char arrays, not const char*
  return vect;
}
/**
 * @brief Allocates a 2D Kokkos view.
 * @tparam T Kokkos view type of rank 2.
 * @param[in] n1 Extent of the first dimension.
 * @param[in] n2 Extent of the second dimension.
 * @return The allocated view.
 */
template <class T>
T allocateArray2D(int n1, int n2) {
  T array(KOKKOSNAME n1, n2);
  return array;
}
/**
 * @brief Allocates a 2D Kokkos view or DualView labeled @p name and, if PRINT_ALLOC_INFO is
 * defined, prints its name and size.
 * @tparam T Kokkos view or DualView type of rank 2.
 * @param[in] n1 Extent of the first dimension.
 * @param[in] n2 Extent of the second dimension.
 * @param[in] name Label of the allocation, see allocateVector().
 * @return The allocated view, zero-initialized.
 */
template <class T>
T allocateArray2D(int n1, int n2, const char *name) {
#ifdef PRINT_ALLOC_INFO
  std::cout << "allocate array : " << name << " of size: (" << n1 << ", " << n2 << ")" << std::endl;
#endif
  T array(std::string(name), n1, n2);
  return array;
}
/**
 * @brief Allocates a 3D Kokkos view.
 * @tparam T Kokkos view type of rank 3.
 * @param[in] n1 Extent of the first dimension.
 * @param[in] n2 Extent of the second dimension.
 * @param[in] n3 Extent of the third dimension.
 * @return The allocated view.
 */
template <class T>
T allocateArray3D(int n1, int n2, int n3) {
  T array(KOKKOSNAME n1, n2, n3);
  return array;
}

/**
 * @brief Wraps an existing device view into a DualView, with a host copy of its content.
 *
 * Both sides are in sync on return. An unallocated view gives an unallocated DualView
 * (extent 0), so "array not provided" stays detectable with extent(0) == 0.
 * @tparam DualViewType Kokkos::DualView type whose t_dev is the type of @p device_view.
 * @param[in] device_view View to wrap; the DualView shares its memory on the device side.
 * @return DualView holding @p device_view and its host copy.
 */
template <class DualViewType>
DualViewType makeDualView(const typename DualViewType::t_dev &device_view) {
  if (!device_view.is_allocated()) return DualViewType();
  auto host_view = Kokkos::create_mirror_view(device_view);
  Kokkos::deep_copy(host_view, device_view);
  return DualViewType(device_view, host_view);
}

/**
 * @brief Debug print of a 3x3 matrix, only for the first two elements (element < 2).
 * @tparam T Type indexable as J[row][col].
 * @param[in] element Element index; nothing is printed if it is 2 or more.
 * @param[in] J 3x3 matrix to print.
 * @param[in] matrixname Name shown in the header line.
 * @param[in] args Values streamed to std::cout before the header line.
 */
template <typename T, typename... Args>
void printJMatrix(const int &element, T &J, std::string matrixname, Args... args) {
  if (element < 2) {
    (std::cout << ... << args) << '\n';
    printf("%s at element %d\n", matrixname.c_str(), element);
    for (int l = 0; l < 3; l++) printf("%f, %f, %f\n", J[l][0], J[l][1], J[l][2]);
  }
}

/**
 * @brief Debug print of a symmetric 3x3 matrix stored in Voigt form, only for the first two
 *        elements (element < 2).
 * @tparam T Type indexable as B[0..5].
 * @param[in] element Element index; nothing is printed if it is 2 or more.
 * @param[in] B Six independent components of the symmetric matrix.
 * @todo VERIFY: is the Voigt order of B (0=xx, 1=yy, 2=zz, 3=yz, 4=xz, 5=xy) the intended one?
 */
template <typename T>
void printBMatrix(const int &element, T &B) {
  if (element < 2) {
    printf("\nB matrix at element %d\n", element);
    printf("%f, %f, %f\n", B[0], B[5], B[4]);
    printf("%f, %f, %f\n", B[5], B[1], B[3]);
    printf("%f, %f, %f\n\n", B[4], B[3], B[2]);
  }
}

/// Declares a std::chrono::system_clock time point named @p timepoint, set to the current time.
#define timewatch(timepoint) \
  std::chrono::time_point<std::chrono::system_clock> timepoint = std::chrono::system_clock::now();
/// Adds the time elapsed since @p starttime to @p accumulatedtime, in clock ticks (period of system_clock, not
/// seconds).
#define accumtime(accumulatedtime, starttime) accumulatedtime += (std::chrono::system_clock::now() - starttime).count();
#endif  // FUNTIDES_CORE_INCLUDE_DATA_TYPE_H_
