#ifndef FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_KIND_H_
#define FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_KIND_H_

namespace solver {
namespace fe {

/**
 * @brief Compile-time selector of the discretization back-end.
 *
 * Used as the key of DiscretizationTraits (Integrals.h), which maps a
 * (polynomial order, kind) pair to the concrete discretization type.
 * @note The solvers do not use it yet; they still select the back-end through
 * IntegralTypeSelector (see docs/design-red-flags.md).
 */
enum class DiscretizationKind {
  kMakutu,        ///< Qk_Hexahedron_Lagrange_GaussLobatto: flat sum-factorization kernels.
  kTensorialGemm  ///< Qk_Hexahedron_Tensorial_GEMM: team-level GEMM kernels.
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_KIND_H_
