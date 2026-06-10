#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
#include "common_macros.h"
namespace solver {
namespace fe {
/**
 * @brief Base Wavefield data structure.
 */
struct Wavefield {
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
  virtual vectorReal getCurrentField(int i) const = 0;

  /**
   * @brief Get the previous field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested current field.
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getPreviousField(int i) const = 0;

  /**
   * @brief Get the previous-previous field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested previous-previous field (empty view if not allocated).
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getPrevPrevField(int i) const = 0;

  /**
   * @brief Check if previous-previous field is allocated.
   * @return True if prevprev buffer exists, false otherwise.
   */
  virtual bool hasPrevPrev() const = 0;

  /**
   * @brief Swap data to advance the wavefield to the next time step.
   * This method exchanges field data pointers.
   * If hasPrevPrev() is true, performs 3-way rotation (prevprev ← prev ← curr ← prevprev).
   * Otherwise performs 2-way swap (curr ↔ prev).
   */
  virtual void swap() = 0;

  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
