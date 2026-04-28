#ifndef FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_

#include <iostream>

#include "data_type.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Data structure for the DG acoustic solver.
 *
 * Holds solution fields (per-element 2-D) and RHS terms (same layout as SEM).
 * Fields are indexed (n_elem, n_dof_per_elem); RHS follows SEM conventions.
 */
struct DGsolverDataAcoustic : public Solver::DataStruct {
  ARRAY_REAL_VIEW pnPrev;     ///< Pressure at previous time step (n_elem, n_dof)
  ARRAY_REAL_VIEW pnCurr;     ///< Pressure at current time step  (n_elem, n_dof)
  ARRAY_REAL_VIEW myRHSTerm;  ///< Source time series (n_rhs, n_sample)
  VECTOR_INT_VIEW rhsElement; ///< Source element indices
  ARRAY_REAL_VIEW rhsWeights; ///< Source weights (n_rhs, n_dof_per_elem)

  bool isDistributed{false};

  DGsolverDataAcoustic(ARRAY_REAL_VIEW prev, ARRAY_REAL_VIEW curr, ARRAY_REAL_VIEW rhsTerm,
                       VECTOR_INT_VIEW rhsElem, ARRAY_REAL_VIEW rhsW)
      : pnPrev(prev), pnCurr(curr), myRHSTerm(rhsTerm), rhsElement(rhsElem), rhsWeights(rhsW) {}

  PROXY_HOST_DEVICE ARRAY_REAL_VIEW getCurrentField(int /*i*/) const { return pnCurr; }
  PROXY_HOST_DEVICE ARRAY_REAL_VIEW getPreviousField(int /*i*/) const { return pnPrev; }
  PROXY_HOST_DEVICE ARRAY_REAL_VIEW getRhsTerm(int /*i*/) const { return myRHSTerm; }
  PROXY_HOST_DEVICE VECTOR_INT_VIEW getRhsElement() const { return rhsElement; }
  PROXY_HOST_DEVICE ARRAY_REAL_VIEW getRhsWeights() const { return rhsWeights; }

  void swapWavefields() { std::swap(pnPrev, pnCurr); }

  void print() const override {
    std::cout << "DGsolverDataAcoustic: " << pnPrev.extent(0) << " elems x " << pnPrev.extent(1) << " dofs"
              << std::endl;
  }
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_
