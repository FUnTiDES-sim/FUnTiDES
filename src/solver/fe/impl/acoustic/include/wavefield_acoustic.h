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

  PROXY_HOST_DEVICE
  WavefieldAcoustic(vectorReal pnGlobalPrev, vectorReal pnGlobalCurr)
      : m_pnGlobalPrev(pnGlobalPrev), m_pnGlobalCurr(pnGlobalCurr) {}

  int getNumFields() const override final { return kNumFields; }

  const char* const* getFieldNames() const override final { return kFieldNames; }

  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const override { return m_pnGlobalCurr; }

  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const override { return m_pnGlobalPrev; }

  void swap() override { std::swap(m_pnGlobalPrev, m_pnGlobalCurr); }

  void swapWithRotation(vectorReal& prevPrevBuffer, int i) override {
    vectorReal tmp = prevPrevBuffer;
    prevPrevBuffer = m_pnGlobalPrev;
    m_pnGlobalPrev = m_pnGlobalCurr;
    m_pnGlobalCurr = tmp;
  }

  void print() const override {
    std::cout << "Pn Global Prev size: " << m_pnGlobalPrev.extent(0) << std::endl;
    std::cout << "Pn Global Curr size: " << m_pnGlobalCurr.extent(0) << std::endl;
  }

  vectorReal m_pnGlobalPrev;  ///< Pressure field at previous time step
  vectorReal m_pnGlobalCurr;  ///< Pressure field at current time step
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
