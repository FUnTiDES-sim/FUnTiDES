#ifndef FUNTIDES_GRADIENT_API_INCLUDE_GRADIENT_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_GRADIENT_H_

#include <string>

#include "data_type.h"

namespace gradient
{

/**
 * @brief Base Gradient data structure.
 */
class Gradient
{
public:
  virtual ~Gradient() = default;

  /**
   * @brief Get the number of gradient fields in this physics.
   * @return The number of fields.
   */
  virtual int getNumGradients() const = 0;

  /**
   * @brief Get the name of a specific gradient field.
   * @param i The index of the gradient field.
   * @return Name of the gradient field.
   */
  virtual std::string getGradientName(int i) const = 0;

  /**
   * @brief Get the field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested field.
   */
  PROXY_HOST_DEVICE
  virtual VECTOR_REAL_VIEW getGradient(int i) const = 0;

  virtual void print() const = 0;
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_API_INCLUDE_GRADIENT_H_
