#ifndef RHS_ACOUSTIC_H_
#define RHS_ACOUSTIC_H_

#include <data_type.h>

#include "rhs.h"

/**
 * @brief Acoustic RHS data structure.
 */
struct RhsAcoustic : public Rhs
{
  RhsAcoustic(ARRAY_REAL_VIEW term, VECTOR_INT_VIEW element,
              ARRAY_REAL_VIEW weights)
      : m_term(term), m_element(element), m_weights(weights)
  {
  }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  ARRAY_REAL_VIEW getTerm(int i) const override { return m_term; }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  VECTOR_INT_VIEW getElement() const { return m_element; }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  ARRAY_REAL_VIEW getWeights() const { return m_weights; }

  void print() const override
  {
    std::cout << "RHS Term size:    " << m_term.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_element.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_weights.extent(0) << std::endl;
  }

  ARRAY_REAL_VIEW m_term;     ///< RHS forcing term
  VECTOR_INT_VIEW m_element;  ///< Source element indices
  ARRAY_REAL_VIEW m_weights;  ///< Forcing weights per node
};

#endif  // RHS_ACOUSTIC_H_