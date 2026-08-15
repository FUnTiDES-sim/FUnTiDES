#ifndef FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_KIND_H_
#define FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_KIND_H_

namespace solver {
namespace fe {

/// Selects a discretization back-end at compile time.
///
/// Mirrors the PhysicsTraits<PHYSICS> pattern: a small enum drives a traits map
/// (DiscretizationTraits, see Integrals.h) instead of threading a bare type
/// template parameter through the whole solver. Zero runtime cost.
enum class DiscretizationKind {
  kMakutu,        ///< Qk_Hexahedron_Lagrange_GaussLobatto  (flat sum-factorization)
  kTensorialGemm  ///< Qk_Hexahedron_Tensorial_GEMM         (team / GEMM)
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_KIND_H_
