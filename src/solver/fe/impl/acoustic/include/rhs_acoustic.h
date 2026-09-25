#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_RHS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_RHS_ACOUSTIC_H_
#include <data_type.h>

#include "rhs.h"

namespace solver {
namespace fe {
/**
 * @brief Source term of the acoustic wave equation (one component).
 *
 * Holds the forcing term, the indices of the elements containing the sources
 * and the interpolation weights of each source on the nodes of its element.
 * Members are Kokkos views: copies are shallow.
 * @todo VERIFY: shapes of m_term, m_element and m_weights (which dimension is
 * the time sample, the source, the node)?
 */
struct RhsAcoustic : public Rhs {
  /// Number of RHS (source) components
  static constexpr int kNumRhsComponents = 1;

  PROXY_HOST_DEVICE RhsAcoustic() = default;
  PROXY_HOST_DEVICE ~RhsAcoustic() = default;
  PROXY_HOST_DEVICE RhsAcoustic(const RhsAcoustic&) = default;
  PROXY_HOST_DEVICE RhsAcoustic& operator=(const RhsAcoustic&) = default;

  /**
   * @brief Builds a source from its views (no copy of the data).
   * @param[in] term Forcing term.
   * @param[in] element Indices of the source elements.
   * @param[in] weights Interpolation weights per node.
   */
  PROXY_HOST_DEVICE
  RhsAcoustic(arrayReal term, vectorInt element, arrayReal weights)
      : m_term(term), m_element(element), m_weights(weights) {}

  /// @return Number of source components (always 1).
  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  /**
   * @brief Returns the forcing term.
   * @param i Component index, ignored (a single component exists).
   */
  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override { return m_term; }

  /// @return Indices of the source elements.
  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_element; }

  /// @return Interpolation weights of the sources.
  PROXY_HOST_DEVICE
  arrayReal getWeights() const { return m_weights; }

  /// @return Interpolation weights of the sources; the component index is ignored.
  PROXY_HOST_DEVICE
  arrayReal getWeights(int) const { return m_weights; }

  /// Prints the extents of the term, element and weights views to stdout.
  void print() const override {
    std::cout << "RHS Term size:    " << m_term.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_element.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_weights.extent(0) << std::endl;
  }

  arrayReal m_term;     ///< Forcing term
  vectorInt m_element;  ///< Indices of the source elements
  arrayReal m_weights;  ///< Interpolation weights per node
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_RHS_ACOUSTIC_H_
