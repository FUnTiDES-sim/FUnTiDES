#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_RHS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_RHS_ACOUSTIC_H_
#include <data_type.h>

#include "rhs.h"

namespace solver {
namespace fe {
/**
 * @brief Acoustic RHS data structure.
 */
struct RhsAcoustic : public Rhs {
  /// Number of RHS (source) components
  static constexpr int kNumRhsComponents = 1;

  // Add explicit device-callable constructors and destructors
  PROXY_HOST_DEVICE RhsAcoustic() = default;
  PROXY_HOST_DEVICE ~RhsAcoustic() = default;
  PROXY_HOST_DEVICE RhsAcoustic(const RhsAcoustic&) = default;
  PROXY_HOST_DEVICE RhsAcoustic& operator=(const RhsAcoustic&) = default;

  PROXY_HOST_DEVICE
  RhsAcoustic(arrayReal term, vectorInt element, arrayReal weights)
      : m_term(term), m_element(element), m_weights(weights) {}

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override { return m_term; }

  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_element; }

  PROXY_HOST_DEVICE
  arrayReal getWeights() const { return m_weights; }

  void print() const override {
    std::cout << "RHS Term size:    " << m_term.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_element.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_weights.extent(0) << std::endl;
  }

  arrayReal m_term;     ///< RHS forcing term
  vectorInt m_element;  ///< Source element indices
  arrayReal m_weights;  ///< Forcing weights per node
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_RHS_ACOUSTIC_H_
