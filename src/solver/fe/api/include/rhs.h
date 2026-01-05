#ifndef RHS_H_
#define RHS_H_

/**
 * @brief Base RHS data structure.
 */
struct Rhs
{
  virtual ~Rhs() = default;

  virtual void print() const = 0;
};

#endif  // RHS_H_