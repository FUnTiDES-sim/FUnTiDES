#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_

#include <string>

#include "data_type.h"

namespace gradient {

/**
 * @brief Read-only set of wavefield arrays (forward or adjoint snapshots) given
 * to a Differentiator, without exposing solver data structures.
 *
 * Each field is a view handle with one value per global node.
 * No caller uses this interface polymorphically; see docs/design-red-flags.md.
 * @see docs/design.md, "Device calls on mesh objects".
 */
class WavefieldView {
 public:
  virtual ~WavefieldView() = default;

  /** @brief Returns the number of fields. */
  virtual int getNumFields() const = 0;

  /**
   * @brief Returns the name of field @p i.
   * @param[in] i Index in [0, getNumFields()).
   */
  virtual std::string getFieldName(int i) const = 0;

  /**
   * @brief Returns the handle of field @p i.
   * @param[in] i Index in [0, getNumFields()).
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getField(int i) const = 0;

  /** @brief Prints a description of the fields to standard output. */
  virtual void print() const = 0;
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_