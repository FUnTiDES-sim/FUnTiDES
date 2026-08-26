#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <typeinfo>
#include <unordered_set>

#include "data_type.h"
#include "dg_penalty.h"
#include "dg_solver_data.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "parallel_topology.h"
#include "physics_traits.h"
#include "physics_traits_acoustic.h"
#include "sem_enums.h"
#include "solver.h"

namespace solver {
namespace fe {

using physicType = utils::enums::physicType;

/**
 * @brief Compile-time face/element DOF lookup tables for a hexahedron of the given order.
 *
 * Depends only on ORDER, not on mesh type, physics, or IS_MODEL_ON_NODES, so it is factored
 * out of DGsolver rather than duplicated per solver instantiation. The indexing math itself
 * lives in model::faceLocalToElemLocal[AtDepth] (face_connectivity.h), which the unstructured
 * mesh builder already calls at runtime and is unit-tested there — this struct just bakes it
 * into a flat array at compile time to avoid the per-lookup switch/call in GPU kernels.
 */
template <int ORDER>
struct DgFaceDofTable {
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  static constexpr int kNumNodesPerFace = (ORDER + 1) * (ORDER + 1);

  static constexpr auto buildFaceToElemDof() {
    std::array<std::array<int, kNumNodesPerFace>, 6> t{};
    for (int f = 0; f < 6; ++f)
      for (int i = 0; i < kNumNodesPerFace; ++i)
        t[f][i] = model::faceLocalToElemLocal(static_cast<model::CubicFace>(f), i, ORDER);
    return t;
  }

  /// @brief Lookup table: kFaceToElemDof[face_id][face_dof_2d] → element-local DOF.
  static constexpr auto kFaceToElemDof = buildFaceToElemDof();

  static constexpr auto buildFaceToElemDofAtDepth() {
    std::array<std::array<std::array<int, ORDER + 1>, kNumNodesPerFace>, 6> t{};
    for (int f = 0; f < 6; ++f)
      for (int i = 0; i < kNumNodesPerFace; ++i)
        for (int m = 0; m <= ORDER; ++m)
          t[f][i][m] = model::faceLocalToElemLocalAtDepth(static_cast<model::CubicFace>(f), i, m, ORDER);
    return t;
  }

  /// @brief Lookup table: kFaceToElemDofAtDepth[face_id][face_dof_2d][depth] → element-local DOF.
  static constexpr auto kFaceToElemDofAtDepth = buildFaceToElemDofAtDepth();
};

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
class DGsolver : public Solver {
 public:
  using Traits = PhysicsTraits<PHYSICS>;

  using DataType = DGsolverDataAcoustic;

  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;

  DGsolver() = default;
  ~DGsolver() = default;

  // --- Implementation of DD Interface ---

  int getNumComponents() const override { return kNumFields; }

  // --- Mandatory overrides for Solver interface ---

  void initFEarrays() override {
    // Here for retrocompatibility
  }

  void allocateFEarrays() override {
    // Here for retrocompatibility
  }

  void initSpongeValues() override {
    // Here for retrocompatibility
  }

  void resetGlobalVectors(int numNodes) override {
    // Here for retrocompatibility, no global vector for DG
  }

  void computeGlobalMassMatrix() override {
    // Here for retrocompatibility, no global mass matrix for DG
  }

  void computeDampingMatrix() override {
    // Here for retrocompatibility, no global damping matrix for DG
  }

  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override {}

  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override;

  vectorReal& getMassMatrixAcoustic() override {
    throw std::runtime_error("getMassMatrixAcoustic not implemented for DG");
  }

  vectorReal& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG");
  }

  vectorReal& getDampingMatrix(int c) override { throw std::runtime_error("getDampingMatrix not implemented for DG"); }

  vectorReal& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
  }

  void setAnisotropyType(model::AnisotropyType type) override {
    // TODO: Implement anisotropy setting
  }

  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    // TODO: Implement SLS attenuation setting
  }

  // --- Core solver methods ---

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override;

  void updateSolutionForward(const float& dt, DataStruct& data) override;
  void updateSolutionBackward(const float& dt, DataStruct& data) override;

  /**
   * @brief Legacy/Serial wrapper.
   * @throws std::runtime_error if called in distributed mode.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override {
    auto& myData = dynamic_cast<DataType&>(data);
    if (myData.isDistributed) {
      throw std::runtime_error(
          "computeOneStep called in distributed mode. Use computeForces() -> "
          "synchronize() -> updateSolutionForward().");
    }
    computeForces(dt, timeSample, data);
    updateSolutionForward(dt, data);
  }

  /**
   * @brief Apply external forcing to the global fields.
   */
  void applyRHSTerm(int timeSample, float dt, const DataType& data);

  /**
   * @brief Update the global solution fields at interior nodes (forward mode).
   */
  void updateFieldsForward(float dt, const DataType& data);

  /**
   * @brief Update the global solution fields (backward/adjoint mode).
   * Note: DG backward mode not yet fully implemented.
   */
  void updateFieldsBackward(float dt, const DataType& data);

  /**
   * @brief Run Verlet update only for a compact subset of elements (forward mode).
   * @param dt        Time step.
   * @param data      Wavefield data.
   * @param elem_list Compact array of element indices to update.
   * @param n_elems   Number of entries in @p elem_list.
   */
  void updateFieldsFromListForward(float dt, const DataType& data, const vectorInt& elem_list, int n_elems);

  /**
   * @brief Run Verlet update only for a compact subset of elements (backward mode).
   * @param dt        Time step.
   * @param data      Wavefield data.
   * @param elem_list Compact array of element indices to update.
   * @param n_elems   Number of entries in @p elem_list.
   * Note: DG backward mode not yet fully implemented.
   */
  void updateFieldsFromListBackward(float dt, const DataType& data, const vectorInt& elem_list, int n_elems);

  /**
   * @brief Kernel 1 — volume mass + SumFact stiffness. Zeros the damping accumulator.
   * @param kNumElem Total number of elements.
   * @param current_field Pressure field at current time step p^n.
   */
  void computeVolumeAndBoundary(int kNumElem, arrayReal current_field);

  /**
   * @brief Kernel 1b — boundary absorbing damping (face-loop, boundary faces only).
   * @param kNumFaces Total number of faces (interior + boundary).
   */
  void computeBoundaryDamping(int kNumFaces);

  /**
   * @brief Kernel 2 — SIPG interface flux terms (reads neighbor fields).
   * @param kNumElem Total number of elements.
   * @param current_field Pressure field at current time step p^n.
   */
  void computeInterfaceFlux(int kNumElem, arrayReal current_field);

  /**
   * @brief Kernel 3 — Verlet time update.
   * @param kNumElem Total number of elements.
   * @param dt Time step size.
   * @param current_field Pressure field at current time step p^n.
   * @param prev_field Pressure field at previous time step p^{n-1}; receives p^{n+1}.
   */
  void applyVerlet(int kNumElem, float dt, arrayReal current_field, arrayReal prev_field);

  /**
   * @brief Given a list of elements, construct the corresponding list of faces for kernels that operate on faces.
   */
  void faceListFromElementList() {
    std::unordered_set<int> visited_faces;
    for (int i = 0; i < m_n_elem_list_; ++i) {
      const int e = m_elem_list_[i];
      for (int f = 0; f < 6; ++f) {
        const int face_id = m_mesh.getGlobalFace(e, static_cast<model::CubicFace>(f));
        visited_faces.insert(face_id);
      }
    }
    m_n_face_list_ = static_cast<int>(visited_faces.size());
    m_face_list_ = allocateVector<vectorInt>(m_n_face_list_, "faceList");
    int i = 0;
    for (const int face_id : visited_faces) {
      m_face_list_(i++) = face_id;
    }
  }

 private:
  MESH_TYPE m_mesh;
  model::FaceConnectivityUnstruct<float, int, ORDER> m_face_connectivity_;
  real_t m_penalty_factor_ = 10.0f;

  // List state used by updateFieldsFromList.
  bool m_list_mode_ = false;
  vectorInt m_elem_list_;
  int m_n_elem_list_ = 0;
  vectorInt m_face_list_;
  int m_n_face_list_ = 0;

  arrayReal m_rhs_elem_;
  arrayReal m_mass_local_;   ///< Per-element mass diagonal (nElem x kPPE)
  arrayReal m_stiff_local_;  ///< Per-element stiffness + interface flux accumulator (nElem x kPPE)
  arrayReal m_damp_local_;   ///< Per-element boundary absorbing damping (nElem x kPPE)

  using DofTable = DgFaceDofTable<ORDER>;
  static constexpr int kPointsPerElement = DofTable::kPointsPerElement;
  static constexpr int knumNodesPerFace = DofTable::kNumNodesPerFace;
  static constexpr auto kFaceToElemDof = DofTable::kFaceToElemDof;
  static constexpr auto kFaceToElemDofAtDepth = DofTable::kFaceToElemDofAtDepth;
};

// Backward Compatibility Aliases
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using DGsolverAcoustic =
    DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kAcoustic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_