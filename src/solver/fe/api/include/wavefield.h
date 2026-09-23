#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
#include "common_macros.h"
namespace solver {
namespace fe {
/**
 * @brief Solution fields of a solver at the time levels of the explicit
 * scheme: current, previous and, in backward mode only, previous-previous.
 *
 * Each field component is a vector with one value per global mesh node.
 *
 * @see docs/design.md, "Time levels and the split time step" for the buffer
 * roles, and "Device calls on mesh objects" for why kernels use the concrete
 * type.
 */
struct Wavefield {
  PROXY_HOST_DEVICE
  virtual ~Wavefield() = default;

  /**
   * @brief Number of field components.
   * @return Bound of the component index of the field getters.
   */
  virtual int getNumFields() const = 0;

  /**
   * @brief Names of the field components, for output and logging.
   * @return Array of getNumFields() null-terminated strings.
   */
  virtual const char* const* getFieldNames() const = 0;

  /**
   * @brief Field at the current time level.
   * @param[in] i Component index, in [0, getNumFields()).
   * @return View on the buffer, one value per global mesh node.
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getCurrentField(int i) const = 0;

  /**
   * @brief Field at the previous time level.
   * @param[in] i Component index, in [0, getNumFields()).
   * @return View on the buffer, one value per global mesh node.
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getPreviousField(int i) const = 0;

  /**
   * @brief Previous-previous buffer, used in backward mode only.
   * @param[in] i Component index, in [0, getNumFields()).
   * @return The buffer, or an empty view when hasPrevPrev() is false.
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getPrevPrevField(int i) const = 0;

  /**
   * @brief Whether the wavefield is in backward mode.
   * @return True if the previous-previous buffer is allocated.
   */
  virtual bool hasPrevPrev() const = 0;

  /**
   * @brief Advance to the next time level by exchanging buffers, without
   * copying data.
   *
   * Call after the solver has written the new time level. Forward mode
   * exchanges current and previous. Backward mode moves previous-previous to
   * current, current to previous and previous to previous-previous.
   */
  virtual void swap() = 0;

  /** @brief Print a debug summary to standard output. */
  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_WAVEFIELD_H_
