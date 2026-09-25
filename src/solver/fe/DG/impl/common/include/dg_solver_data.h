#ifndef FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_

#include <iostream>

#include "data_type.h"
#include "dg_physics_traits_acoustic.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Wavefield and source data handed to the DG acoustic solver at each time step.
 *
 * Bundles the current/previous pressure fields and the right-hand-side source terms.
 * Fields are 2-D arrays indexed (element, dof within element). Both members are
 * lightweight Kokkos view handles, so copies are shallow.
 */
struct DGsolverDataAcoustic : public Solver::DataStruct {
  // Concrete types (not the abstract bases) so that the members can be used on device without virtual dispatch.
  using WavefieldType = typename DGPhysicsTraits::WavefieldType;
  using RhsType = typename DGPhysicsTraits::RhsType;

  WavefieldType m_wavefield;  ///< Current and previous fields, stored by value.
  RhsType m_rhs;              ///< Source terms, stored by value.

  bool isDistributed{false};  ///< @todo VERIFY: meaning of this flag (MPI-distributed run?) and who reads it.

  /**
   * @brief Builds the data from a wavefield and a source.
   * @param[in] wavefield Fields, copied by value (shallow view copy).
   * @param[in] rhs Source terms, copied by value (shallow view copy).
   */
  DGsolverDataAcoustic(const DGWavefieldAcoustic& wavefield, const RhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  /**
   * @brief Returns the current field.
   * @param[in] i Field index.
   * @return View of size (n_elem, n_dof_per_elem).
   */
  PROXY_HOST_DEVICE
  arrayReal getCurrentField(int i) const { return m_wavefield.getCurrentField(i); }

  /**
   * @brief Returns the previous-time-step field.
   * @param[in] i Field index.
   * @return View of size (n_elem, n_dof_per_elem).
   */
  PROXY_HOST_DEVICE
  arrayReal getPreviousField(int i) const { return m_wavefield.getPreviousField(i); }

  /**
   * @brief Returns a source term.
   * @param[in] i Source component index.
   * @return Source term view.
   * @todo VERIFY: shape and unit of the source term view.
   */
  PROXY_HOST_DEVICE
  arrayReal getRhsTerm(int i) const { return m_rhs.getTerm(i); }

  /// @return Indices of the elements that carry a source.
  PROXY_HOST_DEVICE
  vectorInt getRhsElement() const { return m_rhs.getElement(); }

  /// @return Interpolation weights of the sources within their elements.
  PROXY_HOST_DEVICE
  arrayReal getRhsWeights() const { return m_rhs.getWeights(); }

  /// @brief Exchanges the current and previous fields.
  void swapWavefields() { m_wavefield.swap(); }

  /// @brief Prints the extents of the field arrays to stdout.
  void print() const override {
    std::cout << "DGsolverDataAcoustic: " << m_wavefield.getPreviousField(0).extent(0) << " elems x "
              << m_wavefield.getPreviousField(0).extent(1) << " dofs" << std::endl;
  }
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_
