#ifndef SOLVER_FE_API_INCLUDE_GRADIENTS_H_
#define SOLVER_FE_API_INCLUDE_GRADIENTS_H_
namespace solver
{
namespace fe
{
/**
 * @brief Base Gradients data structure.
 */
struct Gradients
{
  virtual ~Gradients() = default;

  /**
   * @brief Get the number of gradient fields in this physics.
   * @return The number of fields.
   */
  virtual int getNumGrads() const = 0;

  /**
   * @brief Get the names of the gradient fields.
   * @return Pointer to array of field names.
   */
  virtual const char* const* getGradsNames() const = 0;

  /**
   * @brief Get the field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested field.
   */
  PROXY_HOST_DEVICE
  virtual VECTOR_REAL_VIEW getCurrentGrads(int i) const = 0;

  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // SOLVER_FE_API_INCLUDE_GRADIENTS_H_
