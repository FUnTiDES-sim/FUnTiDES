#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
#include "common_macros.h"
namespace solver
{
namespace fe
{
/**
 * @brief Base Wavefield data structure.
 */
struct Wavefield
{
  PROXY_HOST_DEVICE
  virtual ~Wavefield() = default;

  /**
   * @brief Get the number of solution fields in this wavefield.
   * @return The number of fields.
   */
  virtual int getNumFields() const = 0;

  /**
   * @brief Get the names of the solution fields.
   * @return Pointer to array of field names.
   */
  virtual const char* const* getFieldNames() const = 0;

  /**
   * @brief Get the current field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested current field.
   */
  PROXY_HOST_DEVICE
  virtual VECTOR_REAL_VIEW getCurrentField(int i) const = 0;

  /**
   * @brief Get the previous field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested current field.
   */
  PROXY_HOST_DEVICE
  virtual VECTOR_REAL_VIEW getPreviousField(int i) const = 0;

  /**
   * @brief Swap data to advance the wavefield to the next time step.
   * This method should exchange the current and previous field data.
   */
  virtual void swap() = 0;

  /**
   * @brief Rotate three buffers without any data copy.
   *
   * Performs a 3-way cyclic rotation of view handles:
   *   prevPrevBuffer ← prev ← curr ← prevPrevBuffer
   *
   * After the call:
   *   - curr      holds the slot previously occupied by prevPrevBuffer
   *               (ready to be overwritten by the next solver step)
   *   - prev      holds what was curr   (the most recently computed field)
   *   - prevPrevBuffer holds what was prev (the field from two steps ago)
   *
   * This is intended for adjoint time-loops where the caller manages one
   * extra external buffer and needs copy-free access to three consecutive
   * time levels for gradient computation.
   *
   * @param prevPrevBuffer  External view handle for the n-2 time level.
   *                        Updated in-place to point to the n-1 level after
   *                        the call.
   * @param i               Field index to swap (0 = first field, 1 = second,
   * etc.).
   */
  virtual void swapWithRotation(VECTOR_REAL_VIEW& prevPrevBuffer, int i) = 0;

  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
