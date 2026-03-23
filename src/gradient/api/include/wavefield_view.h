#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_

namespace gradient {


/**
 * @brief Abstract interface for read-only wavefield views in gradient computation.
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
struct WavefieldView {
  virtual ~WavefieldView() = default;

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
   * @brief Get the field at a specific index.
   * @param i The index of the field to retrieve.
   * @return The requested field.
   */
  PROXY_HOST_DEVICE virtual VECTOR_REAL_VIEW getField(int i) const = 0;

  virtual void print() const = 0;
};

}

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_H_