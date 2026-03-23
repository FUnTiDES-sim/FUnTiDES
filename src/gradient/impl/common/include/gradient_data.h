#ifndef FUNTIDES_GRADIENT_IMPL_COMMON_GRADIENT_DATA_H_
#define FUNTIDES_GRADIENT_IMPL_COMMON_GRADIENT_DATA_H_
#include <iostream>

#include "differentiator.h"
#include "physics_traits_acoustic.h"
#include "physics_traits_elastic.h"

namespace gradient
{

/**
 * @brief Physics-templated container for model parameter gradients.
 *
 * Uses PhysicsTraits to select the concrete
 * gradient type at compile time and stores it by value (lightweight view
 * handles), avoiding virtual dispatch and heap allocation.
 *
 * @tparam PHYSICS The physics type (physicType::kAcoustic or kElastic)
 *
 * Usage:
 *   GradientData<physicType::kAcoustic> data(grad_kappa_view, grad_buoyancy_view);
 *   differentiator->compute(mesh, data, dt);
 *   auto kappa = data.getGradient(0);
 */
template <physicType PHYSICS>
struct GradientData : public Differentiator::DataStruct
{
  using Traits                    = PhysicsTraits<PHYSICS>;

  // Use concrete types from PhysicsTraits to avoid virtual dispatch on device
  using WavefieldViewForwardType  = typename Traits::WavefieldViewForwardType;
  using WavefieldViewBackwardType = typename Traits::WavefieldViewBackwardType;
  using GradientType              = typename Traits::GradientType;

  static constexpr int kNumGradients = GradientType::kNumGrads;

  /**
   * @brief Constructor for acoustic gradient data.
   *
   * @param fwd        Forward wavefield view
   * @param bwd        Adjoint wavefield view
   * @param gradient   Gradient container for acoustic parameters
   */
  template <physicType P = PHYSICS,
            typename    = std::enable_if_t<P == physicType::kAcoustic>>
  GradientData(const WavefieldViewForwardAcoustic&  fwd,
               const WavefieldViewBackwardAcoustic& bwd,
               const GradientAcoustic&          gradient)
      : m_fwd(fwd), m_bwd(bwd), m_gradient(gradient)
  {
  }

  /**
   * @brief Constructor for elastic gradient data.
   *
   * @param fwd        Forward wavefield view
   * @param bwd        Adjoint wavefield view
   * @param gradient   Gradient container for elastic parameters
   */
  template <physicType P = PHYSICS,
            typename    = std::enable_if_t<P == physicType::kElastic>>
  GradientData(const WavefieldViewForwardElastic&  fwd,
               const WavefieldViewBackwardElastic& bwd,
               const GradientElastic&           gradient)
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
    std::cout << "GradientData<" << Traits::kName << ">\n";
    m_fwd.print();
    m_bwd.print();
    m_gradient.print();
  }

  WavefieldViewForwardType  m_fwd;       ///< Forward wavefield snapshot(s)
  WavefieldViewBackwardType m_bwd;       ///< Adjoint wavefield snapshot(s)
  GradientType              m_gradient;  ///< Gradient arrays (view handles)
};

//============================================================================
// Type aliases
//============================================================================

using GradientDataAcoustic = GradientData<physicType::kAcoustic>;
using GradientDataElastic  = GradientData<physicType::kElastic>;

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_COMMON_GRADIENT_DATA_H_
