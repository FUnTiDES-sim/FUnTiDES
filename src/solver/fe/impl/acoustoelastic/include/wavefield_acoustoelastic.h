#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_WAVEFIELD_ACOUSTOELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_WAVEFIELD_ACOUSTOELASTIC_H_

#include <data_type.h>

#include "wavefield.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {

/**
 * @brief Wavefield of a coupled acoustic and elastic domain.
 *
 * Aggregates one acoustic wavefield (pressure) and one elastic wavefield
 * (ux, uy, uz). Fields are addressed by a single index: 0 = pressure,
 * 1 = ux, 2 = uy, 3 = uz.
 */
struct WavefieldAcoustoElastic : public Wavefield {
  /// Total number of solution fields: 1 acoustic (p) + 3 elastic (ux, uy, uz).
  static constexpr int kNumFields = 4;

  /// Field names, in field index order.
  static constexpr const char* kFieldNames[4] = {"pressure", "ux", "uy", "uz"};

  /**
   * @brief Builds a wavefield holding two time levels per field.
   * @param[in] pnGlobalPrev Pressure at the previous time level.
   * @param[in] pnGlobalCurr Pressure at the current time level.
   * @param[in] uxnGlobalPrev x displacement at the previous time level.
   * @param[in] uxnGlobalCurr x displacement at the current time level.
   * @param[in] uynGlobalPrev y displacement at the previous time level.
   * @param[in] uynGlobalCurr y displacement at the current time level.
   * @param[in] uznGlobalPrev z displacement at the previous time level.
   * @param[in] uznGlobalCurr z displacement at the current time level.
   */
  WavefieldAcoustoElastic(vectorReal pnGlobalPrev, vectorReal pnGlobalCurr, vectorReal uxnGlobalPrev,
                          vectorReal uxnGlobalCurr, vectorReal uynGlobalPrev, vectorReal uynGlobalCurr,
                          vectorReal uznGlobalPrev, vectorReal uznGlobalCurr)
      : m_acoustic(pnGlobalPrev, pnGlobalCurr),
        m_elastic(uxnGlobalPrev, uxnGlobalCurr, uynGlobalPrev, uynGlobalCurr, uznGlobalPrev, uznGlobalCurr) {}

  /**
   * @brief Builds a wavefield holding three time levels per field.
   * @param[in] pnGlobalPrevPrev Pressure two time levels back.
   * @param[in] pnGlobalPrev Pressure at the previous time level.
   * @param[in] pnGlobalCurr Pressure at the current time level.
   * @param[in] uxnGlobalPrevPrev x displacement two time levels back.
   * @param[in] uxnGlobalPrev x displacement at the previous time level.
   * @param[in] uxnGlobalCurr x displacement at the current time level.
   * @param[in] uynGlobalPrevPrev y displacement two time levels back.
   * @param[in] uynGlobalPrev y displacement at the previous time level.
   * @param[in] uynGlobalCurr y displacement at the current time level.
   * @param[in] uznGlobalPrevPrev z displacement two time levels back.
   * @param[in] uznGlobalPrev z displacement at the previous time level.
   * @param[in] uznGlobalCurr z displacement at the current time level.
   */
  WavefieldAcoustoElastic(vectorReal pnGlobalPrevPrev, vectorReal pnGlobalPrev, vectorReal pnGlobalCurr,
                          vectorReal uxnGlobalPrevPrev, vectorReal uxnGlobalPrev, vectorReal uxnGlobalCurr,
                          vectorReal uynGlobalPrevPrev, vectorReal uynGlobalPrev, vectorReal uynGlobalCurr,
                          vectorReal uznGlobalPrevPrev, vectorReal uznGlobalPrev, vectorReal uznGlobalCurr)
      : m_acoustic(pnGlobalPrevPrev, pnGlobalPrev, pnGlobalCurr),
        m_elastic(uxnGlobalPrevPrev, uxnGlobalPrev, uxnGlobalCurr, uynGlobalPrevPrev, uynGlobalPrev, uynGlobalCurr,
                  uznGlobalPrevPrev, uznGlobalPrev, uznGlobalCurr) {}

  /// @return Number of solution fields (4).
  int getNumFields() const override final { return kNumFields; }

  /// @return Field names in field index order.
  const char* const* getFieldNames() const override final { return kFieldNames; }

  /**
   * @brief Gets the field at the current time level.
   * @param[in] i Field index: 0 = p, 1 = ux, 2 = uy, 3 = uz.
   * @return The nodal values of the field.
   */
  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const override {
    if (i == 0) return m_acoustic.getCurrentField(0);
    return m_elastic.getCurrentField(i - 1);
  }

  /**
   * @brief Gets the field at the previous time level.
   * @param[in] i Field index: 0 = p, 1 = ux, 2 = uy, 3 = uz.
   * @return The nodal values of the field.
   */
  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const override {
    if (i == 0) return m_acoustic.getPreviousField(0);
    return m_elastic.getPreviousField(i - 1);
  }

  /**
   * @brief Gets the field two time levels back.
   * @param[in] i Field index: 0 = p, 1 = ux, 2 = uy, 3 = uz.
   * @return The nodal values of the field.
   * @pre hasPrevPrev() is true.
   */
  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const override {
    if (i == 0) return m_acoustic.getPrevPrevField(0);
    return m_elastic.getPrevPrevField(i - 1);
  }

  /// @return True if both the acoustic and the elastic parts hold a previous-previous level.
  bool hasPrevPrev() const override {
    return m_acoustic.hasPrevPrev() && m_elastic.hasPrevPrev();
  }

  /// Advances the time levels of both the acoustic and the elastic parts.
  void swap() override {
    m_acoustic.swap();
    m_elastic.swap();
  }

  /// Prints the acoustic part, then the elastic part.
  void print() const override {
    m_acoustic.print();
    m_elastic.print();
  }

  WavefieldAcoustic m_acoustic;  ///< Acoustic pressure wavefield.
  WavefieldElastic m_elastic;    ///< Elastic displacement wavefield.
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_WAVEFIELD_ACOUSTOELASTIC_H_
