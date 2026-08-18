#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
#include <iostream>

#include "differentiator.h"
#include "physics_traits_elastic.h"

namespace gradient {

/**
 * @brief Elastic data passed to a differentiator: forward and adjoint wavefield
 *        views plus the output gradient arrays.
 *
 * The members are lightweight view handles; the class does not own the
 * underlying arrays.
 */
struct DifferentiatorDataElastic : public Differentiator::DataStruct {
  using Traits = PhysicsTraits<utils::enums::physicType::kElastic>;

  using WavefieldViewForwardType = typename Traits::WavefieldViewForwardType;
  using WavefieldViewBackwardType = typename Traits::WavefieldViewBackwardType;
  using GradientType = typename Traits::GradientType;

  /**
   * @brief Builds the data container from the three views.
   *
   * @param fwd           Forward wavefield view
   * @param bwd           Adjoint wavefield view
   * @param gradient      Gradient container for elastic parameters
   * @param firstElement  First element to accumulate (default: all)
   * @param lastElement   One past the last element, or -1 for all
   */
  DifferentiatorDataElastic(const WavefieldViewForwardElastic& fwd, const WavefieldViewBackwardElastic& bwd,
                            const GradientElastic& gradient, int firstElement = 0, int lastElement = -1)
      : m_fwd(fwd), m_bwd(bwd), m_gradient(gradient), m_firstElement(firstElement), m_lastElement(lastElement) {}

  /**
   * @brief Returns the i-th forward wavefield array.
   *
   * @param[in] i  Field index.
   * @return View handle on the forward field.
   * @todo VERIFY: what is the field ordering for index i (which component or derivative does each i select)?
   */
  PROXY_HOST_DEVICE
  vectorReal getForwardField(int i) const { return m_fwd.getField(i); }

  /**
   * @brief Returns the i-th adjoint wavefield array.
   *
   * @param[in] i  Field index.
   * @return View handle on the adjoint field.
   * @todo VERIFY: what is the field ordering for index i (which component or derivative does each i select)?
   */
  PROXY_HOST_DEVICE
  vectorReal getBackwardField(int i) const { return m_bwd.getField(i); }

  /**
   * @brief Returns the i-th gradient array.
   *
   * @param[in] i  Gradient index.
   * @return View handle on the gradient array.
   * @todo VERIFY: which elastic parameter does each index i select (rho, lambda, mu order)?
   */
  PROXY_HOST_DEVICE
  vectorReal getGradient(int i) const { return m_gradient.getGradient(i); }

  /// @brief Prints a description of the three views to stdout. Host only.
  void print() const override {
    std::cout << "DifferentiatorDataElastic\n";
    m_fwd.print();
    m_bwd.print();
    m_gradient.print();
  }

  WavefieldViewForwardType m_fwd;   ///< Forward wavefield snapshot(s)
  WavefieldViewBackwardType m_bwd;  ///< Adjoint wavefield snapshot(s)
  GradientType m_gradient;          ///< Gradient arrays (view handles)
  // Restricts accumulation to a contiguous element range: the coupled
  // acousto-elastic case must gather only from solid elements, and an
  // interface node is shared, so it cannot be masked node-wise afterwards.
  int m_firstElement;
  int m_lastElement;
};

/// Alias of DifferentiatorDataElastic.
using GradientDataElastic = DifferentiatorDataElastic;

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
