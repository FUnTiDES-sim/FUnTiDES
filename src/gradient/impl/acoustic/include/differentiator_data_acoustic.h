#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_DIFFERENTIATOR_DATA_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_DIFFERENTIATOR_DATA_ACOUSTIC_H_
#include <iostream>

#include "differentiator.h"
#include "physics_traits_acoustic.h"

namespace gradient {

/**
 * @brief Data container passed to a differentiator for the acoustic physics.
 *
 * Holds the forward and adjoint wavefield views and the gradient arrays for the
 * acoustic model parameters (kappa, buoyancy). The gradient arrays are the output
 * of the differentiator; the wavefield views are its input.
 *
 * @todo VERIFY: which gradient index is kappa and which is buoyancy, and which field
 * index is which wavefield component in the forward and adjoint views?
 */
struct DifferentiatorDataAcoustic : public Differentiator::DataStruct {
  using Traits = PhysicsTraits<utils::enums::physicType::kAcoustic>;

  using WavefieldViewForwardType = typename Traits::WavefieldViewForwardType;
  using WavefieldViewBackwardType = typename Traits::WavefieldViewBackwardType;
  using GradientType = typename Traits::GradientType;

  /**
   * @brief Constructs the container from the views it stores.
   *
   * @param[in] fwd       Forward wavefield view.
   * @param[in] bwd       Adjoint wavefield view.
   * @param[in] gradient  Gradient container for the acoustic parameters.
   */
  DifferentiatorDataAcoustic(const WavefieldViewForwardAcoustic& fwd, const WavefieldViewBackwardAcoustic& bwd,
                             const GradientAcoustic& gradient)
      : m_fwd(fwd), m_bwd(bwd), m_gradient(gradient) {}

  /// @brief Returns the forward wavefield array number @p i.
  PROXY_HOST_DEVICE
  vectorReal getForwardField(int i) const { return m_fwd.getField(i); }

  /// @brief Returns the adjoint wavefield array number @p i.
  PROXY_HOST_DEVICE
  vectorReal getBackwardField(int i) const { return m_bwd.getField(i); }

  /// @brief Returns the gradient array number @p i.
  PROXY_HOST_DEVICE
  vectorReal getGradient(int i) const { return m_gradient.getGradient(i); }

  /// @brief Prints the forward view, the adjoint view and the gradients to stdout.
  void print() const override {
    std::cout << "DifferentiatorDataAcoustic\n";
    m_fwd.print();
    m_bwd.print();
    m_gradient.print();
  }

  WavefieldViewForwardType m_fwd;   ///< Forward wavefield snapshot(s)
  WavefieldViewBackwardType m_bwd;  ///< Adjoint wavefield snapshot(s)
  GradientType m_gradient;          ///< Gradient arrays (view handles)
};

/// Alias of DifferentiatorDataAcoustic.
using GradientDataAcoustic = DifferentiatorDataAcoustic;

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_DIFFERENTIATOR_DATA_ACOUSTIC_H_
