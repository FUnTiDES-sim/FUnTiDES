#ifndef WAVEFIELD_H_
#define WAVEFIELD_H_

namespace solver
{
namespace fe
{
/**
 * @brief Base Wavefield data structure.
 */
struct Wavefield
{
  virtual ~Wavefield() = default;

  /**
   * @brief Get the current field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested current field.
   */
#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  virtual VECTOR_REAL_VIEW getCurrentField(int i) const = 0;

  /**
   * @brief Get the previous field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested current field.
   */
#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  virtual VECTOR_REAL_VIEW getPreviousField(int i) const = 0;

  /**
   * @brief Swap data to advance the wavefield to the next time step.
   * This method should exchange the current and previous field data.
   */
  virtual void swap() = 0;

  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver

#endif  // WAVEFIELD_H_