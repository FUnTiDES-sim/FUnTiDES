#ifndef RHS_H_
#define RHS_H_

/**
 * @brief Base RHS data structure.
 */
struct Rhs
{
  virtual ~Rhs() = default;

  /**
   * @brief Get a specific RHS term.
   * @param i The index of the term to retrieve.
   * @return The requested RHS term.
   */
  virtual ARRAY_REAL_VIEW getTerm(int i) const = 0;

  /**
   * @brief Get the element indices associated with the RHS.
   * @return The element indices.
   */
  virtual VECTOR_INT_VIEW getElement() const = 0;

  /**
   * @brief Get the weights associated with the RHS.
   * @return The weights.
   */
  virtual ARRAY_REAL_VIEW getWeights() const = 0;

  virtual void print() const = 0;
};

#endif  // RHS_H_