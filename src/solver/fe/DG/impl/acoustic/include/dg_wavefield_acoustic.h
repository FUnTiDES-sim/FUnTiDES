#ifndef FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_WAVEFIELD_ACOUSTIC_H_
#include <data_type.h>

namespace solver {
namespace fe {
/**
 * @brief Acoustic wavefield data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct DGWavefieldAcoustic {
  /// Number of solution fields (1 for pressure)
  static constexpr int kNumFields = 1;

  /// Primary field name
  static constexpr const char* kFieldNames[1] = {"pressure"};

  // Add explicit device-callable constructors and destructors
  PROXY_HOST_DEVICE DGWavefieldAcoustic() = default;
  PROXY_HOST_DEVICE ~DGWavefieldAcoustic() = default;
  PROXY_HOST_DEVICE DGWavefieldAcoustic(const DGWavefieldAcoustic&) = default;
  PROXY_HOST_DEVICE DGWavefieldAcoustic& operator=(const DGWavefieldAcoustic&) = default;

  PROXY_HOST_DEVICE
  DGWavefieldAcoustic(arrayReal pnPrev, arrayReal pnCurr) : m_pnPrev(pnPrev), m_pnCurr(pnCurr) {}

  int getNumFields() const { return kNumFields; }

  const char* const* getFieldNames() const { return kFieldNames; }

  PROXY_HOST_DEVICE
  arrayReal getCurrentField(int i) const { return m_pnCurr; }

  PROXY_HOST_DEVICE
  arrayReal getPreviousField(int i) const { return m_pnPrev; }

  void swap() { std::swap(m_pnPrev, m_pnCurr); }

  void rotate(arrayReal& prevPrevBuffer, int i) {
    arrayReal tmp = prevPrevBuffer;
    prevPrevBuffer = m_pnPrev;
    m_pnPrev = m_pnCurr;
    m_pnCurr = tmp;
  }

  void print() const {
    std::cout << "Pn Prev size: " << m_pnPrev.extent(0) << " elems " << m_pnPrev.extent(1) << " dofs" << std::endl;
    std::cout << "Pn Curr size: " << m_pnCurr.extent(0) << " elems " << m_pnCurr.extent(1) << " dofs" << std::endl;
  }

  arrayReal m_pnPrev;  ///< Pressure field at previous time step (n_elem, n_dof)
  arrayReal m_pnCurr;  ///< Pressure field at current time step (n_elem, n_dof)
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_WAVEFIELD_ACOUSTIC_H_
