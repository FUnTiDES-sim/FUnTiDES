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
  arrayReal pnPrev;      ///< Pressure at previous time step (n_elem, n_dof)
  arrayReal pnCurr;      ///< Pressure at current time step  (n_elem, n_dof)
  arrayReal myRHSTerm;   ///< Source time series (n_rhs, n_sample)
  vectorInt rhsElement;  ///< Source element indices
  arrayReal rhsWeights;  ///< Source weights (n_rhs, n_dof_per_elem)

  bool isDistributed{false};

  DGsolverDataAcoustic(arrayReal prev, arrayReal curr, arrayReal rhsTerm, vectorInt rhsElem, arrayReal rhsW)
      : pnPrev(prev), pnCurr(curr), myRHSTerm(rhsTerm), rhsElement(rhsElem), rhsWeights(rhsW) {}

  PROXY_HOST_DEVICE arrayReal getCurrentField(int /*i*/) const { return pnCurr; }
  PROXY_HOST_DEVICE arrayReal getPreviousField(int /*i*/) const { return pnPrev; }
  PROXY_HOST_DEVICE arrayReal getRhsTerm(int /*i*/) const { return myRHSTerm; }
  PROXY_HOST_DEVICE vectorInt getRhsElement() const { return rhsElement; }
  PROXY_HOST_DEVICE arrayReal getRhsWeights() const { return rhsWeights; }

  void swapWavefields() { std::swap(pnPrev, pnCurr); }

  void print() const override {
    std::cout << "DGsolverDataAcoustic: " << pnPrev.extent(0) << " elems x " << pnPrev.extent(1) << " dofs"
              << std::endl;
  }
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_
