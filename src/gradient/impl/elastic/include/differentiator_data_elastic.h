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
   * @param[in] fwd       Forward wavefield view.
   * @param[in] bwd       Adjoint wavefield view.
   * @param[in] gradient  Gradient container for the elastic parameters.
   */
  DifferentiatorDataElastic(const WavefieldViewForwardElastic& fwd, const WavefieldViewBackwardElastic& bwd,
                            const GradientElastic& gradient)
      : m_fwd(fwd), m_bwd(bwd), m_gradient(gradient) {}

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

  WavefieldViewForwardType m_fwd;   ///< Forward wavefield view.
  WavefieldViewBackwardType m_bwd;  ///< Adjoint wavefield view.
  GradientType m_gradient;          ///< Gradient arrays (view handles).
};

/// Alias of DifferentiatorDataElastic.
using GradientDataElastic = DifferentiatorDataElastic;

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
