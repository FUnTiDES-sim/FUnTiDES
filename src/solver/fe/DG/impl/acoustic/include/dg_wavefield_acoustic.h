#ifndef FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_WAVEFIELD_ACOUSTIC_H_
#include <data_type.h>

namespace solver {
namespace fe {
/**
 * @brief Pressure wavefield of the acoustic discontinuous Galerkin solver.
 *
 * Holds the pressure at the previous and current time steps as flat 2D arrays
 * indexed by (element, local dof).
 */
struct DGWavefieldAcoustic {
  static constexpr int kNumFields = 1;  ///< Number of solution fields (pressure only).

  static constexpr const char* kFieldNames[1] = {"pressure"};  ///< Name of each field.

  PROXY_HOST_DEVICE DGWavefieldAcoustic() = default;
  PROXY_HOST_DEVICE ~DGWavefieldAcoustic() = default;
  PROXY_HOST_DEVICE DGWavefieldAcoustic(const DGWavefieldAcoustic&) = default;
  PROXY_HOST_DEVICE DGWavefieldAcoustic& operator=(const DGWavefieldAcoustic&) = default;

  /**
   * @brief Wraps existing pressure arrays without copying (shallow view copy).
   * @param[in] pnPrev Pressure at the previous time step, shape (n_elem, n_dof).
   * @param[in] pnCurr Pressure at the current time step, shape (n_elem, n_dof).
   */
  PROXY_HOST_DEVICE
  DGWavefieldAcoustic(arrayReal pnPrev, arrayReal pnCurr) : m_pnPrev(pnPrev), m_pnCurr(pnCurr) {}

  /// @return Number of solution fields.
  int getNumFields() const { return kNumFields; }

  /// @return Array of kNumFields field names.
  const char* const* getFieldNames() const { return kFieldNames; }

  /// @return Pressure at the current time step. The index is ignored.
  PROXY_HOST_DEVICE
  arrayReal getCurrentField(int i) const { return m_pnCurr; }

  /// @return Pressure at the previous time step. The index is ignored.
  PROXY_HOST_DEVICE
  arrayReal getPreviousField(int i) const { return m_pnPrev; }

  /// @brief Exchanges the previous and current pressure arrays.
  void swap() { std::swap(m_pnPrev, m_pnCurr); }

  /// @brief Prints the extents of both arrays to stdout.
  void print() const {
    std::cout << "Pn Prev size: " << m_pnPrev.extent(0) << " elems " << m_pnPrev.extent(1) << " dofs" << std::endl;
    std::cout << "Pn Curr size: " << m_pnCurr.extent(0) << " elems " << m_pnCurr.extent(1) << " dofs" << std::endl;
  }

  arrayReal m_pnPrev;  ///< Pressure at the previous time step, shape (n_elem, n_dof).
  arrayReal m_pnCurr;  ///< Pressure at the current time step, shape (n_elem, n_dof).
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_WAVEFIELD_ACOUSTIC_H_
