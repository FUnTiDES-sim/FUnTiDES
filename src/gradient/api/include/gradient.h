#ifndef FUNTIDES_GRADIENT_API_INCLUDE_GRADIENT_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_GRADIENT_H_

#include <string>

#include "data_type.h"

namespace gradient {

/**
 * @brief Set of gradient arrays of one physics, one array per model
 * parameter, held as view handles.
 *
 * Each array holds one value per global node or per element, depending on
 * Differentiator::isModelOnNodes().
 * No caller uses this interface polymorphically; see docs/design-red-flags.md.
 * @see docs/design.md, "Device calls on mesh objects".
 */
class Gradient {
 public:
  virtual ~Gradient() = default;

  /** @brief Returns the number of gradient arrays. */
  virtual int getNumGradients() const = 0;

  /**
   * @brief Returns the name of gradient array @p i.
   * @param[in] i Index in [0, getNumGradients()).
   */
  virtual std::string getGradientName(int i) const = 0;

  /**
   * @brief Returns the handle of gradient array @p i.
   * @param[in] i Index in [0, getNumGradients()).
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getGradient(int i) const = 0;

  /** @brief Prints a description of the arrays to standard output. */
  virtual void print() const = 0;
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_API_INCLUDE_GRADIENT_H_
