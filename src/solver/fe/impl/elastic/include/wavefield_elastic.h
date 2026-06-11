#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_WAVEFIELD_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_WAVEFIELD_ELASTIC_H_
#include <data_type.h>

#include "wavefield.h"

namespace solver {
namespace fe {
/**
 * @brief Elastic wavefield data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct WavefieldElastic : public Wavefield {
  /// Field names for each component
  static constexpr const char* kFieldNames[3] = {"ux", "uy", "uz"};
  /// Number of solution fields (3 for displacement vector)
  static constexpr int kNumFields = 3;

  // Add explicit device-callable constructors and destructors
  PROXY_HOST_DEVICE WavefieldElastic() = default;
  PROXY_HOST_DEVICE ~WavefieldElastic() = default;
  PROXY_HOST_DEVICE WavefieldElastic(const WavefieldElastic&) = default;
  PROXY_HOST_DEVICE WavefieldElastic& operator=(const WavefieldElastic&) = default;

  /*
   *  @brief Constructor for forward simulation (2-buffer mode).
   *  Contains current and previous fields for each displacement component.
   */
  PROXY_HOST_DEVICE
  WavefieldElastic(vectorReal uxnGlobalPrev, vectorReal uxnGlobalCurr, vectorReal uynGlobalPrev,
                   vectorReal uynGlobalCurr, vectorReal uznGlobalPrev, vectorReal uznGlobalCurr)
      : m_uxnGlobalPrevPrev(),
        m_uxnGlobalPrev(uxnGlobalPrev),
        m_uxnGlobalCurr(uxnGlobalCurr),
        m_uynGlobalPrevPrev(),
        m_uynGlobalPrev(uynGlobalPrev),
        m_uynGlobalCurr(uynGlobalCurr),
        m_uznGlobalPrevPrev(),
        m_uznGlobalPrev(uznGlobalPrev),
        m_uznGlobalCurr(uznGlobalCurr) {}

  /*
   *  @brief Constructor for adjoint/backward simulation (3-buffer mode).
   *  Contains current, previous, and previous-previous fields for each displacement component.
   */
  PROXY_HOST_DEVICE
  WavefieldElastic(vectorReal uxnGlobalPrevPrev, vectorReal uxnGlobalPrev, vectorReal uxnGlobalCurr,
                   vectorReal uynGlobalPrevPrev, vectorReal uynGlobalPrev, vectorReal uynGlobalCurr,
                   vectorReal uznGlobalPrevPrev, vectorReal uznGlobalPrev, vectorReal uznGlobalCurr)
      : m_uxnGlobalPrevPrev(uxnGlobalPrevPrev),
        m_uxnGlobalPrev(uxnGlobalPrev),
        m_uxnGlobalCurr(uxnGlobalCurr),
        m_uynGlobalPrevPrev(uynGlobalPrevPrev),
        m_uynGlobalPrev(uynGlobalPrev),
        m_uynGlobalCurr(uynGlobalCurr),
        m_uznGlobalPrevPrev(uznGlobalPrevPrev),
        m_uznGlobalPrev(uznGlobalPrev),
        m_uznGlobalCurr(uznGlobalCurr) {}

  int getNumFields() const override final { return kNumFields; }

  const char* const* getFieldNames() const override final { return kFieldNames; }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const override {
    switch (i) {
      case 0:
        return m_uxnGlobalCurr;
      case 1:
        return m_uynGlobalCurr;
      case 2:
        return m_uznGlobalCurr;
      default:
        return m_uxnGlobalCurr;  // make it cuda happy
    }
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const override {
    switch (i) {
      case 0:
        return m_uxnGlobalPrev;
      case 1:
        return m_uynGlobalPrev;
      case 2:
        return m_uznGlobalPrev;
      default:
        return m_uxnGlobalCurr;  // make it cuda happy
    }
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const override {
    switch (i) {
      case 0:
        return m_uxnGlobalPrevPrev;
      case 1:
        return m_uynGlobalPrevPrev;
      case 2:
        return m_uznGlobalPrevPrev;
      default:
        return m_uxnGlobalPrevPrev;  // make it cuda happy
    }
  }

  bool hasPrevPrev() const override { return m_uxnGlobalPrevPrev.extent(0) > 0; }

  void swap() override {
    if (hasPrevPrev()) {
      // Backward mode: curr ← prevPrev (new value), prev ← curr, prevPrev ← prev
      vectorReal tempUx = m_uxnGlobalCurr;
      vectorReal tempUy = m_uynGlobalCurr;
      vectorReal tempUz = m_uznGlobalCurr;

      m_uxnGlobalCurr = m_uxnGlobalPrevPrev;  // New value goes to curr
      m_uynGlobalCurr = m_uynGlobalPrevPrev;
      m_uznGlobalCurr = m_uznGlobalPrevPrev;

      m_uxnGlobalPrevPrev = m_uxnGlobalPrev;
      m_uynGlobalPrevPrev = m_uynGlobalPrev;
      m_uznGlobalPrevPrev = m_uznGlobalPrev;

      m_uxnGlobalPrev = tempUx;
      m_uynGlobalPrev = tempUy;
      m_uznGlobalPrev = tempUz;
    } else {
      // 2-way swap: curr ↔ prev
      std::swap(m_uxnGlobalPrev, m_uxnGlobalCurr);
      std::swap(m_uynGlobalPrev, m_uynGlobalCurr);
      std::swap(m_uznGlobalPrev, m_uznGlobalCurr);
    }
  }

  void print() const override {
    std::cout << "Ux Global Prev size: " << m_uxnGlobalPrev.extent(0) << std::endl;
    std::cout << "Ux Global Curr size: " << m_uxnGlobalCurr.extent(0) << std::endl;
    std::cout << "Uy Global Prev size: " << m_uynGlobalPrev.extent(0) << std::endl;
    std::cout << "Uy Global Curr size: " << m_uynGlobalCurr.extent(0) << std::endl;
    std::cout << "Uz Global Prev size: " << m_uznGlobalPrev.extent(0) << std::endl;
    std::cout << "Uz Global Curr size: " << m_uznGlobalCurr.extent(0) << std::endl;
    if (hasPrevPrev()) {
      std::cout << "Ux Global PrevPrev size: " << m_uxnGlobalPrevPrev.extent(0) << std::endl;
      std::cout << "Uy Global PrevPrev size: " << m_uynGlobalPrevPrev.extent(0) << std::endl;
      std::cout << "Uz Global PrevPrev size: " << m_uznGlobalPrevPrev.extent(0) << std::endl;
    }
  }

  vectorReal m_uxnGlobalPrevPrev;  ///< Displacement field in x at n-2 time step (for adjoint)
  vectorReal m_uxnGlobalPrev;      ///< Displacement field in x at previous time step
  vectorReal m_uxnGlobalCurr;      ///< Displacement field in x at current time step
  vectorReal m_uynGlobalPrevPrev;  ///< Displacement field in y at n-2 time step (for adjoint)
  vectorReal m_uynGlobalPrev;      ///< Displacement field in y at previous time step
  vectorReal m_uynGlobalCurr;      ///< Displacement field in y at current time step
  vectorReal m_uznGlobalPrevPrev;  ///< Displacement field in z at n-2 time step (for adjoint)
  vectorReal m_uznGlobalPrev;      ///< Displacement field in z at previous time step
  vectorReal m_uznGlobalCurr;      ///< Displacement field in z at current time step
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_WAVEFIELD_ELASTIC_H_
