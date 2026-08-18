#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
#include <iostream>

#include "differentiator.h"
#include "physics_traits_elastic.h"

namespace gradient {

/**
 * @brief Elastic data container for differentiator computation.
 *
 * Stores forward and adjoint wavefield views along with gradient arrays
 * for elastic model parameters (rho, lambda, mu). Passed to
 * DifferentiatorElastic::compute() at runtime.
 *
 * Usage:
 *   DifferentiatorDataElastic data(fwd, bwd, gradient);
 *   float dt = 0.001;
 *   differentiator->compute(mesh, data, dt);
 *   auto rho = data.getGradient(0);
 */
struct DifferentiatorDataElastic : public Differentiator::DataStruct {
  using Traits = PhysicsTraits<utils::enums::physicType::kElastic>;

  using WavefieldViewForwardType = typename Traits::WavefieldViewForwardType;
  using WavefieldViewBackwardType = typename Traits::WavefieldViewBackwardType;
  using GradientType = typename Traits::GradientType;

  /**
   * @brief Construct elastic differentiator data.
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

  PROXY_HOST_DEVICE
  vectorReal getForwardField(int i) const { return m_fwd.getField(i); }

  PROXY_HOST_DEVICE
  vectorReal getBackwardField(int i) const { return m_bwd.getField(i); }

  PROXY_HOST_DEVICE
  vectorReal getGradient(int i) const { return m_gradient.getGradient(i); }

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

using GradientDataElastic = DifferentiatorDataElastic;

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_DIFFERENTIATOR_DATA_ELASTIC_H_
