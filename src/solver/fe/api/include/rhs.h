#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_RHS_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_RHS_H_
#include "common_macros.h"
namespace solver
{
namespace fe
{
/**
 * @brief Base RHS data structure.
 */
struct Rhs
{
  PROXY_HOST_DEVICE
  virtual ~Rhs() = default;

  /**
   * @brief Get the number of RHS components.
   * @return The number of RHS components.
   */
  virtual int getNumRhsComponents() const = 0;

  /**
   * @brief Get a specific RHS term.
   * @param i The index of the term to retrieve.
   * @return The requested RHS term.
   */
  PROXY_HOST_DEVICE
  virtual ARRAY_REAL_VIEW getTerm(int i) const = 0;

  /**
   * @brief Get the element indices associated with the RHS.
   * @return The element indices.
   */
  PROXY_HOST_DEVICE
  virtual VECTOR_INT_VIEW getElement() const = 0;

  /**
   * @brief Get the weights associated with the RHS.
   * @return The weights.
   */
  PROXY_HOST_DEVICE
  virtual ARRAY_REAL_VIEW getWeights() const = 0;

  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_RHS_H_
