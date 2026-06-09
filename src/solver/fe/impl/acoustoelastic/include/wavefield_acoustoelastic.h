#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_WAVEFIELD_ACOUSTOELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_WAVEFIELD_ACOUSTOELASTIC_H_

#include <data_type.h>

#include "wavefield.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {

/**
 * @brief Combined wavefield for the acousto-elastic coupled solver.
 *
 * Holds both the acoustic pressure field and the three elastic displacement
 * components at two consecutive time levels (previous and current).
 * This struct is passed to SEMsolverAcoustoElastic at each time step.
 */
struct WavefieldAcoustoElastic : public Wavefield {
  /// Total number of solution fields: 1 acoustic (p) + 3 elastic (ux, uy, uz)
  static constexpr int kNumFields = 4;

  /// Field names in order: p, ux, uy, uz
  static constexpr const char* kFieldNames[4] = {"pressure", "ux", "uy", "uz"};

  WavefieldAcoustoElastic(vectorReal pnGlobalPrev, vectorReal pnGlobalCurr, vectorReal uxnGlobalPrev,
                          vectorReal uxnGlobalCurr, vectorReal uynGlobalPrev, vectorReal uynGlobalCurr,
                          vectorReal uznGlobalPrev, vectorReal uznGlobalCurr)
      : m_acoustic(pnGlobalPrev, pnGlobalCurr),
        m_elastic(uxnGlobalPrev, uxnGlobalCurr, uynGlobalPrev, uynGlobalCurr, uznGlobalPrev, uznGlobalCurr) {}

  int getNumFields() const override final { return kNumFields; }

  const char* const* getFieldNames() const override final { return kFieldNames; }

  /**
   * @brief Get the current field by index.
   *
   * Index mapping: 0=p, 1=ux, 2=uy, 3=uz.
   */
  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const override {
    if (i == 0) return m_acoustic.getCurrentField(0);
    return m_elastic.getCurrentField(i - 1);
  }

  /**
   * @brief Get the previous field by index.
   *
   * Index mapping: 0=p, 1=ux, 2=uy, 3=uz.
   */
  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const override {
    if (i == 0) return m_acoustic.getPreviousField(0);
    return m_elastic.getPreviousField(i - 1);
  }

  void swap() override {
    m_acoustic.swap();
    m_elastic.swap();
  }

  void rotate(vectorReal& prevPrevBuffer, int i) override {
    if (i == 0) {
      m_acoustic.rotate(prevPrevBuffer, 0);
      return;
    }
    m_elastic.rotate(prevPrevBuffer, i - 1);
  }

  void print() const override {
    m_acoustic.print();
    m_elastic.print();
  }

  WavefieldAcoustic m_acoustic;  ///< Acoustic pressure wavefield
  WavefieldElastic m_elastic;    ///< Elastic displacement wavefield
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_WAVEFIELD_ACOUSTOELASTIC_H_
