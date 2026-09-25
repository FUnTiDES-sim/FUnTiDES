/**
 * @file
 * @brief Compile-time selection of the hexahedral discretization type from a
 * polynomial order and a back-end.
 */
#ifndef FUNTIDES_DISCRETIZATION_FE_IMPL_COMMON_INTEGRALS_H_
#define FUNTIDES_DISCRETIZATION_FE_IMPL_COMMON_INTEGRALS_H_

#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"
#include "Qk_Hexahedron_Tensorial.h"
#include "fe_discretization_kind.h"

namespace solver {
namespace fe {

/**
 * @brief Maps a polynomial order and a back-end to the discretization type.
 *
 * Each specialization provides `type`, the discretization class, and
 * `kHasGemm`, true when `type` provides the team-level GEMM kernels, so that a
 * caller can branch with `if constexpr`.
 * @tparam ORDER Polynomial order, from 1 to 9.
 * @tparam KIND Discretization back-end.
 * @note Currently unused by the solvers (see docs/design-red-flags.md).
 */
template <int ORDER, DiscretizationKind KIND>
struct DiscretizationTraits;

template <int ORDER>
struct DiscretizationTraits<ORDER, DiscretizationKind::kMakutu> {
  using type = typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
  static constexpr bool kHasGemm = false;
};

template <int ORDER>
struct DiscretizationTraits<ORDER, DiscretizationKind::kTensorialGemm> {
  using type = typename Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<ORDER>::type;
  static constexpr bool kHasGemm = true;
};

/**
 * @brief Discretization type for polynomial order ORDER and back-end KIND.
 * @tparam ORDER Polynomial order, from 1 to 9.
 * @tparam KIND Discretization back-end.
 */
template <int ORDER, DiscretizationKind KIND>
using DiscretizationType = typename DiscretizationTraits<ORDER, KIND>::type;

}  // namespace fe
}  // namespace solver

/**
 * @brief Maps a polynomial order and an IntegralType value to the
 * discretization type, exposed as `type`.
 * @tparam ORDER Polynomial order, from 1 to 9.
 * @tparam METHOD_TYPE One of the IntegralType values.
 * @deprecated Use solver::fe::DiscretizationTraits (see
 * docs/design-red-flags.md).
 */
template <int ORDER, int METHOD_TYPE>
struct IntegralTypeSelector;

/**
 * @brief Back-end keys for the METHOD_TYPE parameter of IntegralTypeSelector.
 */
namespace IntegralType {
enum { MAKUTU, TENSORIAL_GEMM };
}

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::MAKUTU> {
  using type = typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
};

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::TENSORIAL_GEMM> {
  using type = typename Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<ORDER>::type;
};

#endif  // FUNTIDES_DISCRETIZATION_FE_IMPL_COMMON_INTEGRALS_H_
