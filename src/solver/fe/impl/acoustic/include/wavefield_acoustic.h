#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
#include <data_type.h>

#include "wavefield.h"

namespace solver {
namespace fe {
/**
 * @brief Pressure wavefield holding up to three time levels of one scalar field.
 *
 * The buffers are flat vectors. The n-2 buffer is optional: it is empty when
 * the object is built with two buffers.
 * @todo VERIFY: size and indexing of the pressure vectors (one value per global node?).
 */
struct WavefieldAcoustic : public Wavefield {
  /// Number of solution fields (1: pressure).
  static constexpr int kNumFields = 1;

  /// Names of the solution fields, indexed like the fields.
  static constexpr const char* kFieldNames[1] = {"pressure"};

  PROXY_HOST_DEVICE WavefieldAcoustic() = default;
  PROXY_HOST_DEVICE ~WavefieldAcoustic() = default;
  PROXY_HOST_DEVICE WavefieldAcoustic(const WavefieldAcoustic&) = default;
  PROXY_HOST_DEVICE WavefieldAcoustic& operator=(const WavefieldAcoustic&) = default;

  /**
   * @brief Builds a wavefield with three time levels (n-2, n-1, n).
   * @param[in] pnGlobalPrevPrev Pressure at time step n-2.
   * @param[in] pnGlobalPrev Pressure at time step n-1.
   * @param[in] pnGlobalCurr Pressure at time step n.
   */
  PROXY_HOST_DEVICE
  WavefieldAcoustic(vectorReal pnGlobalPrevPrev, vectorReal pnGlobalPrev, vectorReal pnGlobalCurr)
      : m_pnGlobalPrevPrev(pnGlobalPrevPrev), m_pnGlobalPrev(pnGlobalPrev), m_pnGlobalCurr(pnGlobalCurr) {}

  /**
   * @brief Builds a wavefield with two time levels (n-1, n); the n-2 buffer stays empty.
   * @param[in] pnGlobalPrev Pressure at time step n-1.
   * @param[in] pnGlobalCurr Pressure at time step n.
   */
  PROXY_HOST_DEVICE
  WavefieldAcoustic(vectorReal pnGlobalPrev, vectorReal pnGlobalCurr)
      : m_pnGlobalPrevPrev(), m_pnGlobalPrev(pnGlobalPrev), m_pnGlobalCurr(pnGlobalCurr) {}

  /// @brief Returns the number of solution fields.
  int getNumFields() const override final { return kNumFields; }

  /// @brief Returns the field names, indexed like the fields.
  const char* const* getFieldNames() const override final { return kFieldNames; }

  /**
   * @brief Returns the pressure at the current time step.
   * @param[in] i Field index, ignored (there is a single field).
   */
  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const override { return m_pnGlobalCurr; }

  /**
   * @brief Returns the pressure at the previous time step.
   * @param[in] i Field index, ignored (there is a single field).
   */
  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const override { return m_pnGlobalPrev; }

  /**
   * @brief Returns the pressure at time step n-2 (empty if unused).
   * @param[in] i Field index, ignored (there is a single field).
   */
  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const override { return m_pnGlobalPrevPrev; }

  /// @brief Returns true if the n-2 buffer is allocated.
  bool hasPrevPrev() const override { return m_pnGlobalPrevPrev.extent(0) > 0; }

  /**
   * @brief Rotates the time levels after a step, without copying data.
   *
   * With three buffers, the n-2 buffer becomes current (it is expected to hold
   * the newly computed values), current becomes previous and previous becomes n-2.
   * With two buffers, current and previous are exchanged.
   */
  void swap() override {
    if (hasPrevPrev()) {
      vectorReal temp = m_pnGlobalCurr;
      m_pnGlobalCurr = m_pnGlobalPrevPrev;
      m_pnGlobalPrevPrev = m_pnGlobalPrev;
      m_pnGlobalPrev = temp;
    } else {
      std::swap(m_pnGlobalPrev, m_pnGlobalCurr);
    }
  }

  /// @brief Prints the size of each allocated buffer to stdout.
  void print() const override {
    std::cout << "Pn Global Prev size: " << m_pnGlobalPrev.extent(0) << std::endl;
    std::cout << "Pn Global Curr size: " << m_pnGlobalCurr.extent(0) << std::endl;
    if (hasPrevPrev()) {
      std::cout << "Pn Global PrevPrev size: " << m_pnGlobalPrevPrev.extent(0) << std::endl;
    }
  }

  vectorReal m_pnGlobalPrevPrev;  ///< Pressure at time step n-2; empty when unused.
  vectorReal m_pnGlobalPrev;      ///< Pressure at time step n-1.
  vectorReal m_pnGlobalCurr;      ///< Pressure at time step n.
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
