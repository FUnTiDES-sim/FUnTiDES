#ifndef FUNTIDES_DISCRETIZATION_FE_IMPL_COMMON_INTEGRALS_H_
#define FUNTIDES_DISCRETIZATION_FE_IMPL_COMMON_INTEGRALS_H_

#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"
#include "Qk_Hexahedron_Tensorial.h"
#include "fe_discretization_kind.h"

namespace solver {
namespace fe {

// Maps a (polynomial order, back-end) pair to the concrete discretization type,
// mirroring the PhysicsTraits<PHYSICS> policy pattern. `kHasGemm` exposes the
// team/GEMM capability as plain data so the solver can dispatch with
// `if constexpr` instead of a SFINAE probe.
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

template <int ORDER, DiscretizationKind KIND>
using DiscretizationType = typename DiscretizationTraits<ORDER, KIND>::type;

}  // namespace fe
}  // namespace solver

// ---------------------------------------------------------------------------
// Legacy selection API (DEPRECATED). Kept so the still-INTEGRAL_TYPE-templated
// SEMsolver keeps compiling during the migration to DiscretizationKind.
// Delete once the solver signature is switched over.
// ---------------------------------------------------------------------------
template <int ORDER, int METHOD_TYPE>
struct IntegralTypeSelector;

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
