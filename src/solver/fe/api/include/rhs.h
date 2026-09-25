#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_RHS_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_RHS_H_
#include "common_macros.h"
namespace solver {
namespace fe {
/**
 * @brief Point sources (right-hand side) injected by a solver at each time
 * step: per source, its element, its time function per component and its
 * weights on the element DOFs.
 *
 * @see docs/design.md, "Source terms" for the array layouts, and "Device calls
 * on mesh objects" for why kernels use the concrete type.
 */
struct Rhs {
  PROXY_HOST_DEVICE
  virtual ~Rhs() = default;

  /**
   * @brief Number of source components.
   * @return Bound of the component index of getTerm().
   */
  virtual int getNumRhsComponents() const = 0;

  /**
   * @brief Time function of one component, for every source.
   * @param[in] i Component index, in [0, getNumRhsComponents()).
   * @return Array of shape (number of sources, number of time samples).
   */
  PROXY_HOST_DEVICE
  virtual arrayReal getTerm(int i) const = 0;

  /**
   * @brief Element containing each source.
   * @return One element index per source.
   */
  PROXY_HOST_DEVICE
  virtual vectorInt getElement() const = 0;

  /**
   * @brief Weights that spread each source over the DOFs of its element.
   * @return Array of shape (number of sources, DOFs per element), indexed by
   * element-local DOF. Only one weight set is exposed for all components.
   */
  PROXY_HOST_DEVICE
  virtual arrayReal getWeights() const = 0;

  /** @brief Print a debug summary to standard output. */
  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_RHS_H_
