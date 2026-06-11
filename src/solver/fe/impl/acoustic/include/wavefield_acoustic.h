#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
#include <data_type.h>

#include "wavefield.h"

namespace solver {
namespace fe {
/**
 * @brief Acoustic wavefield data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct WavefieldAcoustic : public Wavefield {
  /// Number of solution fields (1 for pressure)
  static constexpr int kNumFields = 1;

  /// Primary field name
  static constexpr const char* kFieldNames[1] = {"pressure"};

  // Add explicit device-callable constructors and destructors
  PROXY_HOST_DEVICE WavefieldAcoustic() = default;
  PROXY_HOST_DEVICE ~WavefieldAcoustic() = default;
  PROXY_HOST_DEVICE WavefieldAcoustic(const WavefieldAcoustic&) = default;
  PROXY_HOST_DEVICE WavefieldAcoustic& operator=(const WavefieldAcoustic&) = default;

  /*
   *  @brief Constructor for forward simulation (2-buffer mode).
   *  Contains current and previous fields for each displacement component.
   */
  PROXY_HOST_DEVICE
  WavefieldAcoustic(vectorReal pnGlobalPrevPrev, vectorReal pnGlobalPrev, vectorReal pnGlobalCurr)
      : m_pnGlobalPrevPrev(pnGlobalPrevPrev), m_pnGlobalPrev(pnGlobalPrev), m_pnGlobalCurr(pnGlobalCurr) {}

  /*
   *  @brief Constructor for adjoint/backward simulation (3-buffer mode).
   *  Contains current, previous, and previous-previous fields for each displacement component.
   */
  PROXY_HOST_DEVICE
  WavefieldAcoustic(vectorReal pnGlobalPrev, vectorReal pnGlobalCurr)
      : m_pnGlobalPrevPrev(), m_pnGlobalPrev(pnGlobalPrev), m_pnGlobalCurr(pnGlobalCurr) {}

  int getNumFields() const override final { return kNumFields; }

  const char* const* getFieldNames() const override final { return kFieldNames; }

  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const override { return m_pnGlobalCurr; }

  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const override { return m_pnGlobalPrev; }

  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const override { return m_pnGlobalPrevPrev; }

  bool hasPrevPrev() const override { return m_pnGlobalPrevPrev.extent(0) > 0; }

  void swap() override {
    if (hasPrevPrev()) {
      // Backward mode: curr ← prevPrev (new value), prev ← curr, prevPrev ← prev
      vectorReal temp = m_pnGlobalCurr;
      m_pnGlobalCurr = m_pnGlobalPrevPrev;  // New value goes to curr
      m_pnGlobalPrevPrev = m_pnGlobalPrev;
      m_pnGlobalPrev = temp;
    } else {
      // 2-way swap: curr ↔ prev
      std::swap(m_pnGlobalPrev, m_pnGlobalCurr);
    }
  }

  void print() const override {
    std::cout << "Pn Global Prev size: " << m_pnGlobalPrev.extent(0) << std::endl;
    std::cout << "Pn Global Curr size: " << m_pnGlobalCurr.extent(0) << std::endl;
    if (hasPrevPrev()) {
      std::cout << "Pn Global PrevPrev size: " << m_pnGlobalPrevPrev.extent(0) << std::endl;
    }
  }

  vectorReal m_pnGlobalPrevPrev;  ///< Pressure field at n-2 time step (for adjoint)
  vectorReal m_pnGlobalPrev;      ///< Pressure field at previous time step
  vectorReal m_pnGlobalCurr;      ///< Pressure field at current time step
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
