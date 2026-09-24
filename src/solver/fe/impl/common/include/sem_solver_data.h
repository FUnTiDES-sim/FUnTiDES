#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_DATA_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_DATA_H_
#include <iostream>

#include "physics_traits_acoustic.h"
#include "physics_traits_elastic.h"
#include "solver.h"

namespace solver {
namespace fe {

using physicType = utils::enums::physicType;

/**
 * @brief Solution fields and source terms of a SEM solver, for one physics.
 *
 * Holds the wavefield and the source (RHS) by value, as lightweight view
 * handles that can be copied into device kernels. The number of fields and of
 * RHS components is fixed at compile time by the physics.
 *
 * @tparam PHYSICS Physics type (kAcoustic or kElastic).
 */
template <physicType PHYSICS>
struct SEMsolverData : public Solver::DataStruct {
  using Traits = PhysicsTraits<PHYSICS>;
  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;  ///< Number of wavefield components.
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;    ///< Number of RHS components.

  using WavefieldType = typename Traits::WavefieldType;  ///< Concrete wavefield type (no virtual dispatch on device).
  using RhsType = typename Traits::RhsType;              ///< Concrete source type (no virtual dispatch on device).

  /**
   * @brief Builds the data of an acoustic solver (enabled for kAcoustic only).
   * @param[in] wavefield Acoustic wavefield.
   * @param[in] rhs Acoustic source.
   */
  template <physicType P = PHYSICS, typename = std::enable_if_t<P == physicType::kAcoustic>>
  SEMsolverData(const WavefieldAcoustic& wavefield, const RhsAcoustic& rhs) : m_wavefield(wavefield), m_rhs(rhs) {}

  /**
   * @brief Builds the data of an elastic solver (enabled for kElastic only).
   * @param[in] wavefield Elastic wavefield.
   * @param[in] rhs Elastic source.
   * @param[in] isDistributed Value stored in the isDistributed flag.
   */
  template <physicType P = PHYSICS, typename = std::enable_if_t<P == physicType::kElastic>>
  SEMsolverData(const WavefieldElastic& wavefield, const RhsElastic& rhs, const bool isDistributed = false)
      : m_wavefield(wavefield), m_rhs(rhs), isDistributed(isDistributed) {}

  /**
   * @brief Source term of one RHS component.
   * @param[in] i RHS component index.
   * @return View of the source term.
   */
  PROXY_HOST_DEVICE
  arrayReal getRhsTerm(int i) const { return m_rhs.getTerm(i); }

  /**
   * @brief Element indices of the sources.
   * @return View of the source element indices.
   */
  PROXY_HOST_DEVICE
  vectorInt getRhsElement() const { return m_rhs.getElement(); }

  /**
   * @brief Interpolation weights of the sources, all components together.
   * @return View of the source weights.
   */
  PROXY_HOST_DEVICE
  arrayReal getRhsWeights() const { return m_rhs.getWeights(); }

  /**
   * @brief Interpolation weights of the sources for one component.
   * @param[in] i Component index.
   * @return View of the source weights.
   */
  PROXY_HOST_DEVICE
  arrayReal getRhsWeights(int i) const { return m_rhs.getWeights(i); }

  /**
   * @brief Wavefield component at the current time step.
   * @param[in] i Component index.
   * @return View of size numNodes.
   * @todo VERIFY: is the view size numNodes for every physics?
   */
  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const { return m_wavefield.getCurrentField(i); }

  /**
   * @brief Wavefield component at the previous time step.
   * @param[in] i Component index.
   * @return View of the field.
   */
  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const { return m_wavefield.getPreviousField(i); }

  /**
   * @brief Wavefield component two time steps back.
   * @param[in] i Component index.
   * @return View of the field.
   */
  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const { return m_wavefield.getPrevPrevField(i); }

  /// @brief Rotates the time-level buffers of the wavefield.
  void swapWavefields() { m_wavefield.swap(); }

  /// @brief Prints the physics name and the sizes of fields and sources to stdout.
  void print() const override {
    std::cout << "SEMsolverData<" << Traits::kName << ">" << std::endl;
    for (int f = 0; f < kNumFields; ++f) {
      std::cout << "Field[" << f << "] (" << Traits::WavefieldType::kFieldNames[f]
                << ") size: " << getCurrentField(f).extent(0) << std::endl;
    }
    for (int r = 0; r < kNumRhs; ++r) {
      std::cout << "RHS[" << r << "] size: " << getRhsTerm(r).extent(0) << std::endl;
    }
    std::cout << "RHS Element size: " << getRhsElement().extent(0) << std::endl;
    std::cout << "RHS Weights size: " << getRhsWeights().extent(0) << std::endl;
  }

  bool isDistributed{false};  ///< @todo VERIFY: meaning of this flag (set by the elastic constructor only).
  WavefieldType m_wavefield;  ///< Wavefield, stored by value (view handles).
  RhsType m_rhs;              ///< Source, stored by value (view handles).
};

using SEMsolverDataAcoustic = SEMsolverData<physicType::kAcoustic>;  ///< Data of the acoustic SEM solver.
using SEMsolverDataElastic = SEMsolverData<physicType::kElastic>;    ///< Data of the elastic SEM solver.

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_DATA_H_
