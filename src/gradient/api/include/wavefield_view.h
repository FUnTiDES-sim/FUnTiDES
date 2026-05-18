#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_

#include <string>

#include "data_type.h"

namespace gradient {

/**
 * @brief Abstract interface for read-only wavefield views (aka snapshots) in
 * gradient computation.
 *
 * WavefieldView provides a physics-agnostic interface to access the required
 * forward and adjoint wavefield snapshots for gradient computation, without
 * exposing solver-specific data structures.
 *
 * Concrete implementations (e.g. WavefieldViewForwardAcoustic) will wrap solver
 * wavefield data and expose only the fields needed by the gradient kernels.
 *
 * This design allows the Differentiator to operate on abstract wavefield views,
 * decoupling it from solver internals and enabling flexible data management.
 */
class WavefieldView {
 public:
  virtual ~WavefieldView() = default;

  /**
   * @brief Get the number of solution fields in this wavefield.
   * @return The number of fields.
   */
  virtual int getNumFields() const = 0;

  /**
   * @brief Get the name of a specific solution field.
   * @param i The index of the field.
   * @return Name of the field.
   */
  virtual std::string getFieldName(int i) const = 0;

  /**
   * @brief Get the field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested field.
   */
  PROXY_HOST_DEVICE
  virtual vectorReal getField(int i) const = 0;

  virtual void print() const = 0;
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_