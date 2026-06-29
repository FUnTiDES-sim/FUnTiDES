#ifndef FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_H_
#define FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_H_

#include <type_traits>

namespace discretization {
namespace fe {
namespace api {

/// Empty tag base for every Qk discretization back-end.
///
/// Lets the solver and the tests constrain on "is a discretization" WITHOUT a
/// vtable: all real dispatch stays compile-time, which is mandatory on GPU.
struct FeDiscretizationTag {};

/// Compile-time contract ("named requirements") for a discretization type.
///
/// C++17 has no concepts, so the contract is expressed as static_asserts.
/// Call it from a unit test or a static_assert site to document and enforce the
/// static surface every back-end must provide. The heavy static kernels
/// (jacobianTransformation, computeBMatrix, computeMassTerm, computeDampingTerm,
/// computeStiffnessTermSumFact, ...) are exercised by the shared TYPED_TEST
/// suite -- a single set of tests for all back-ends.
template <typename T>
constexpr bool AssertFeDiscretization() {
  static_assert(std::is_base_of<FeDiscretizationTag, T>::value,
                "discretization must derive from FeDiscretizationTag");
  static_assert(T::num1dNodes > 0, "missing/invalid num1dNodes");
  static_assert(T::numNodes == T::num1dNodes * T::num1dNodes * T::num1dNodes,
                "numNodes must equal num1dNodes^3");
  static_assert(T::numQuadraturePoints == T::numNodes,
                "GLL: number of quadrature points must equal number of nodes");
  static_assert(T::maxSupportPoints == T::numNodes, "maxSupportPoints mismatch");
  return true;
}

}  // namespace api
}  // namespace fe
}  // namespace discretization
#endif  // FUNTIDES_DISCRETIZATION_FE_API_FE_DISCRETIZATION_H_
