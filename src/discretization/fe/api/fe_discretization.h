#ifndef FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_H_
#define FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_H_

#include <type_traits>

namespace discretization {
namespace fe {
namespace api {

/**
 * @brief Empty base that marks a type as a Qk hexahedral discretization
 * back-end.
 *
 * Carries no virtual function: back-ends are selected at compile time because
 * virtual dispatch is unusable in device code.
 * @see docs/design.md, "Device calls on mesh objects".
 * @note No back-end currently derives from this tag (see
 * docs/design-red-flags.md).
 */
struct FeDiscretizationTag {};

/**
 * @brief Checks at compile time the static members every discretization
 * back-end must provide.
 *
 * Compilation fails unless @p T derives from FeDiscretizationTag, has
 * num1dNodes > 0, numNodes == num1dNodes^3, and numQuadraturePoints and
 * maxSupportPoints both equal to numNodes (Gauss-Lobatto quadrature on the
 * support points). The static kernels are not checked.
 *
 * @tparam T Candidate discretization type.
 * @return Always true, so the call can sit inside a static_assert.
 * @note Not called anywhere yet (see docs/design-red-flags.md).
 */
template <typename T>
constexpr bool AssertFeDiscretization() {
  static_assert(std::is_base_of<FeDiscretizationTag, T>::value, "discretization must derive from FeDiscretizationTag");
  static_assert(T::num1dNodes > 0, "missing/invalid num1dNodes");
  static_assert(T::numNodes == T::num1dNodes * T::num1dNodes * T::num1dNodes, "numNodes must equal num1dNodes^3");
  static_assert(T::numQuadraturePoints == T::numNodes, "GLL: number of quadrature points must equal number of nodes");
  static_assert(T::maxSupportPoints == T::numNodes, "maxSupportPoints mismatch");
  return true;
}

}  // namespace api
}  // namespace fe
}  // namespace discretization
#endif  // FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_H_
