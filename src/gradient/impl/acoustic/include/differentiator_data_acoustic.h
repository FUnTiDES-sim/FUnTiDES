#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_DIFFERENTIATOR_DATA_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_DIFFERENTIATOR_DATA_ACOUSTIC_H_
#include <iostream>

#include "differentiator.h"
#include "physics_traits_acoustic.h"

namespace gradient
{

/**
 * @brief Acoustic data container for differentiator computation.
 *
 * Stores forward and adjoint wavefield views along with gradient arrays
 * for acoustic model parameters (kappa, buoyancy). Passed to
 * DifferentiatorAcoustic::compute() at runtime.
 *
 * Usage:
 *   DifferentiatorDataAcoustic data(fwd, bwd, gradient);
 *   differentiator->compute(mesh, data);
 *   auto kappa = data.getGradient(0);
 */
struct DifferentiatorDataAcoustic : public Differentiator::DataStruct
{
  using Traits = PhysicsTraits<utils::enums::physicType::kAcoustic>;

  using WavefieldViewForwardType = typename Traits::WavefieldViewForwardType;
  using WavefieldViewBackwardType = typename Traits::WavefieldViewBackwardType;
  using GradientType = typename Traits::GradientType;

  /**
   * @brief Construct acoustic differentiator data.
   *
   * @param fwd       Forward wavefield view
   * @param bwd       Adjoint wavefield view
   * @param gradient  Gradient container for acoustic parameters
   */
  DifferentiatorDataAcoustic(const WavefieldViewForwardAcoustic& fwd,
                             const WavefieldViewBackwardAcoustic& bwd,
                             const GradientAcoustic& gradient)
      : m_fwd(fwd), m_bwd(bwd), m_gradient(gradient)
  {
  }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getForwardField(int i) const { return m_fwd.getField(i); }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getBackwardField(int i) const { return m_bwd.getField(i); }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getGradient(int i) const
  {
    return m_gradient.getGradient(i);
  }

  void print() const override
  {
    std::cout << "DifferentiatorDataAcoustic\n";
    m_fwd.print();
    m_bwd.print();
    m_gradient.print();
  }

  WavefieldViewForwardType m_fwd;   ///< Forward wavefield snapshot(s)
  WavefieldViewBackwardType m_bwd;  ///< Adjoint wavefield snapshot(s)
  GradientType m_gradient;          ///< Gradient arrays (view handles)
};

using GradientDataAcoustic = DifferentiatorDataAcoustic;

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_DIFFERENTIATOR_DATA_ACOUSTIC_H_
