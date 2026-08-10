#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_RHS_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_RHS_ELASTIC_H_
#include <data_type.h>

#include "rhs.h"

namespace solver {
namespace fe {
/**
 * @brief Elastic RHS data structure.
 */
struct RhsElastic : public Rhs {
  /// Number of RHS (source) components
  static constexpr int kNumRhsComponents = 3;

  /// Default constructor: produces an empty (no-forcing) RHS.
  /// All views are null/zero-extent; safe as long as RHS is not accessed.
  RhsElastic() = default;

  RhsElastic(arrayReal termx, arrayReal termy, arrayReal termz, vectorInt element, arrayReal weights)
      : m_termx(termx),
        m_termy(termy),
        m_termz(termz),
        m_element(element),
        m_weights(weights),
        m_weightsy(weights),
        m_weightsz(weights) {}

  /// Variant with one weight array per component, for a source whose spatial
  /// pattern differs between components (moment tensor, explosion in a solid
  /// expressed as a pressure gradient).
  RhsElastic(arrayReal termx, arrayReal termy, arrayReal termz, vectorInt element, arrayReal weightsx,
             arrayReal weightsy, arrayReal weightsz)
      : m_termx(termx),
        m_termy(termy),
        m_termz(termz),
        m_element(element),
        m_weights(weightsx),
        m_weightsy(weightsy),
        m_weightsz(weightsz) {}

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override {
    switch (i) {
      case 0:
        return m_termx;
      case 1:
        return m_termy;
      case 2:
        return m_termz;
      default:
        return m_termx;  // make it cuda happy
    }
  }

  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_element; }

  PROXY_HOST_DEVICE
  arrayReal getWeights() const { return m_weights; }

  PROXY_HOST_DEVICE
  arrayReal getWeights(int i) const {
    switch (i) {
      case 1:
        return m_weightsy;
      case 2:
        return m_weightsz;
      default:
        return m_weights;
    }
  }

  void print() const override {
    std::cout << "RHSx Term size:   " << m_termx.extent(0) << std::endl;
    std::cout << "RHSy Term size:   " << m_termy.extent(0) << std::endl;
    std::cout << "RHSz Term size:   " << m_termz.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_element.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_weights.extent(0) << std::endl;
  }

  arrayReal m_termx;    ///< X-component forcing term
  arrayReal m_termy;    ///< Y-component forcing term
  arrayReal m_termz;    ///< Z-component forcing term
  vectorInt m_element;  ///< Source element indices
  arrayReal m_weights;  ///< Forcing weights per node, X component
  arrayReal m_weightsy;  ///< Forcing weights per node, Y component
  arrayReal m_weightsz;  ///< Forcing weights per node, Z component
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_RHS_ELASTIC_H_