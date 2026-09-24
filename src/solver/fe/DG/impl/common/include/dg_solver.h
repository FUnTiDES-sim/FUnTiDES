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
 * @brief Compile-time face/element DOF lookup tables for a hexahedron of a given order.
 *
 * Depends only on ORDER, not on mesh type, physics or model storage. The tables bake
 * model::faceLocalToElemLocal and model::faceLocalToElemLocalAtDepth (face_connectivity.h)
 * into flat arrays so that GPU kernels avoid a per-lookup call.
 *
 * @tparam ORDER Polynomial order of the element.
 */
template <int ORDER>
struct DgFaceDofTable {
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);  ///< DOFs per element.
  static constexpr int kNumNodesPerFace = (ORDER + 1) * (ORDER + 1);                 ///< DOFs per face.

  /// @brief Builds kFaceToElemDof at compile time.
  static constexpr auto buildFaceToElemDof() {
    std::array<std::array<int, kNumNodesPerFace>, 6> t{};
    for (int f = 0; f < 6; ++f)
      for (int i = 0; i < kNumNodesPerFace; ++i)
        t[f][i] = model::faceLocalToElemLocal(static_cast<model::CubicFace>(f), i, ORDER);
    return t;
  }

  /// Element-local DOF of a face DOF: kFaceToElemDof[face_id][face_dof_2d]. face_id follows model::CubicFace.
  static constexpr auto kFaceToElemDof = buildFaceToElemDof();

  /// @brief Builds kFaceToElemDofAtDepth at compile time.
  static constexpr auto buildFaceToElemDofAtDepth() {
    std::array<std::array<std::array<int, ORDER + 1>, kNumNodesPerFace>, 6> t{};
    for (int f = 0; f < 6; ++f)
      for (int i = 0; i < kNumNodesPerFace; ++i)
        for (int m = 0; m <= ORDER; ++m)
          t[f][i][m] = model::faceLocalToElemLocalAtDepth(static_cast<model::CubicFace>(f), i, m, ORDER);
    return t;
  }

  /// Element-local DOF at a given depth below a face DOF: kFaceToElemDofAtDepth[face_id][face_dof_2d][depth].
  /// @todo VERIFY: depth 0 is the face itself and depth ORDER the opposite face?
  static constexpr auto kFaceToElemDofAtDepth = buildFaceToElemDofAtDepth();
};

/**
 * @brief Discontinuous Galerkin solver on unstructured hexahedral meshes (SIPG fluxes, Verlet time stepping).
 *
 * Fields are stored per element (no global assembly), so the global-array hooks of Solver are
 * no-ops. A time step is split into computeForces() (element and face kernels) and
 * updateSolutionForward() (Verlet update), so that a caller can synchronize in between.
 *
 * @tparam ORDER Polynomial order.
 * @tparam INTEGRAL_TYPE Integral back-end.
 * @tparam MESH_TYPE Concrete mesh type.
 * @tparam IS_MODEL_ON_NODES True if model properties are stored per node, false if per element.
 * @tparam PHYSICS Physics of the solver.
 */
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

  /// @brief Number of wavefield components.
  int getNumComponents() const override { return kNumFields; }

  /// @brief No-op: DG has no global arrays to initialize.
  void initFEarrays() override {}

  /// @brief No-op: DG has no global arrays to allocate.
  void allocateFEarrays() override {}

  /// @brief No-op: sponge values are not used by DG.
  void initSpongeValues() override {}

  /// @brief No-op: DG has no global vectors.
  void resetGlobalVectors(int numNodes) override {}

  /// @brief No-op: DG has no global mass matrix.
  void computeGlobalMassMatrix() override {}

  /// @brief No-op: DG has no global damping matrix.
  void computeDampingMatrix() override {}

  /// @brief No-op: output of a vector field is not supported.
  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override {}

  /// @brief Outputs the per-element field @p field at time sample @p t under the name @p fieldName.
  /// @todo VERIFY: meaning of the in,out argument @p e.
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override;

  /// @throws std::runtime_error always: there is no global mass matrix.
  vectorReal& getMassMatrixAcoustic() override {
    throw std::runtime_error("getMassMatrixAcoustic not implemented for DG");
  }

  /// @throws std::runtime_error always: there is no global mass matrix.
  vectorReal& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG");
  }

  /// @throws std::runtime_error always: there is no global damping matrix.
  vectorReal& getDampingMatrix(int c) override { throw std::runtime_error("getDampingMatrix not implemented for DG"); }

  /// @throws std::runtime_error always: there is no global force vector.
  vectorReal& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
  }

  /// @brief Ignored: anisotropy is not implemented.
  void setAnisotropyType(model::AnisotropyType type) override {}

  /// @brief Ignored: attenuation is not implemented.
  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {}

  /**
   * @brief Builds the per-element arrays and face connectivity from the mesh.
   * @param[in] mesh Mesh, expected to be of type MESH_TYPE.
   * @param[in] sponge_size Ignored.
   * @param[in] surface_sponge Ignored.
   * @param[in] taper_delta Ignored.
   */
  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  /**
   * @brief Computes the volume, boundary damping and interface flux terms of one time step.
   * @param[in] dt Time step.
   * @param[in] timeSample Index of the current time sample, used to evaluate the source.
   * @param[in,out] data Must be a DataType.
   */
  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override;

  /// @brief Advances the fields one step with the terms of computeForces().
  void updateSolutionForward(const float& dt, DataStruct& data) override;

  /// @brief Backward (adjoint) update.
  /// @note Backward mode is not fully implemented for DG.
  void updateSolutionBackward(const float& dt, DataStruct& data) override;

  /**
   * @brief Runs a full time step: computeForces() followed by updateSolutionForward().
   * @throws std::runtime_error if @p data is flagged as distributed (the caller must then
   * synchronize between the two calls).
   * @throws std::bad_cast if @p data is not a DataType.
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
   * @brief Adds the source term to the global fields.
   * @param timeSample Index of the current time sample.
   * @param dt Time step.
   * @param data Wavefield and source data.
   */
  void applyRHSTerm(int timeSample, float dt, const DataType& data);

  /**
   * @brief Verlet update of the global fields at interior nodes (forward mode).
   * @param dt Time step.
   * @param data Wavefield data.
   */
  void updateFieldsForward(float dt, const DataType& data);

  /**
   * @brief Verlet update of the global fields (backward mode).
   * @param dt Time step.
   * @param data Wavefield data.
   * @note Backward mode is not fully implemented for DG.
   */
  void updateFieldsBackward(float dt, const DataType& data);

  /**
   * @brief Verlet update restricted to a subset of elements (forward mode).
   * @param dt Time step.
   * @param data Wavefield data.
   * @param elem_list Compact array of element indices to update.
   * @param n_elems Number of entries in @p elem_list.
   */
  void updateFieldsFromListForward(float dt, const DataType& data, const vectorInt& elem_list, int n_elems);

  /**
   * @brief Verlet update restricted to a subset of elements (backward mode).
   * @param dt Time step.
   * @param data Wavefield data.
   * @param elem_list Compact array of element indices to update.
   * @param n_elems Number of entries in @p elem_list.
   * @note Backward mode is not fully implemented for DG.
   */
  void updateFieldsFromListBackward(float dt, const DataType& data, const vectorInt& elem_list, int n_elems);

  /**
   * @brief Kernel 1: volume mass and sum-factorized stiffness terms. Zeroes the damping accumulator.
   * @param kNumElem Total number of elements.
   * @param current_field Pressure field at the current time step p^n.
   */
  void computeVolumeAndBoundary(int kNumElem, arrayReal current_field);

  /**
   * @brief Kernel 2: boundary absorbing damping and SIPG interface flux terms, fused in one face loop.
   *
   * The two terms are mutually exclusive per face and use disjoint accumulators.
   *
   * @param kNumFaces Total number of faces (interior and boundary).
   * @param current_field Pressure field at the current time step p^n.
   * @note One thread per face. A team-per-face variant was measured on GH200 at order 6: it is
   * 3x faster on a 20^3 mesh (25k faces) but slower on 40^3 (197k faces) and 100x45x60
   * (823k faces). It only wins when one thread per face cannot fill the device.
   */
  void computeBoundaryDampingAndInterfaceFlux(int kNumFaces, arrayReal current_field);

  /**
   * @brief Kernel 3: Verlet time update.
   * @param kNumElem Total number of elements.
   * @param dt Time step.
   * @param current_field Pressure field at the current time step p^n.
   * @param prev_field Pressure field at the previous time step p^{n-1}; receives p^{n+1}.
   */
  void applyVerlet(int kNumElem, float dt, arrayReal current_field, arrayReal prev_field);

  /**
   * @brief Fills m_face_list_ and m_n_face_list_ with the faces of the elements in m_elem_list_.
   *
   * Each face appears once. The order of the resulting list is unspecified.
   */
  void faceListFromElementList() {
    std::unordered_set<int> visited_faces;
    for (int i = 0; i < m_n_elem_list_; ++i) {
      const int e = m_elem_list_[i];
      for (int f = 0; f < 6; ++f) {
        const int face_id = m_face_connectivity_.getGlobalFace(e, static_cast<model::CubicFace>(f));
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

  /// @brief SIPG penalty factor.
  real_t getPenaltyFactor() const { return m_penalty_factor_; }

  /**
   * @brief Injects an externally built face connectivity, so that face ids are shared with the caller.
   *
   * Call before computeFEInit(), which then skips its own build().
   *
   * @param face_connectivity Already-built connectivity for the same mesh.
   */
  void setFaceConnectivity(const model::FaceConnectivityUnstruct<float, int, ORDER>& face_connectivity) {
    m_face_connectivity_ = face_connectivity;
  }

  /**
   * @brief Rebuilds the face connectivity sampling face DOFs at the mesh geometric order instead of ORDER.
   *
   * Needed when ORDER is lower than the order of the mesh: otherwise the faces on the plus side
   * (kXPlus, kYPlus, kZPlus) place the face-normal coordinate at ORDER instead of the true far
   * edge, and the neighbor node matching fails silently on those faces. No-op in effect when
   * the mesh order equals ORDER.
   *
   * @param geometricOrder Polynomial order of the mesh (mesh.getOrder()).
   */
  void rebuildFaceConnectivityGeometry(int geometricOrder) { m_face_connectivity_.build(m_mesh, geometricOrder); }

 private:
  MESH_TYPE m_mesh;
  model::FaceConnectivityUnstruct<float, int, ORDER> m_face_connectivity_;
  real_t m_penalty_factor_ = 12.0f;

 public:
  bool m_list_mode_ = false;  ///< If true, the update kernels run on m_elem_list_ and m_face_list_ only.
  vectorInt m_elem_list_;     ///< Element ids of the subset, first m_n_elem_list_ entries valid.
  int m_n_elem_list_ = 0;     ///< Number of valid entries of m_elem_list_.
  vectorInt m_face_list_;     ///< Face ids of the subset, first m_n_face_list_ entries valid.
  int m_n_face_list_ = 0;     ///< Number of valid entries of m_face_list_.

  arrayReal m_rhs_elem_;     ///< Per-element source term.
  arrayReal m_mass_local_;   ///< Per-element mass diagonal (nElem x kPPE)
  arrayReal m_stiff_local_;  ///< Per-element stiffness + interface flux accumulator (nElem x kPPE)
  arrayReal m_damp_local_;   ///< Per-element boundary absorbing damping (nElem x kPPE)

  using DofTable = DgFaceDofTable<ORDER>;
  static constexpr int kPointsPerElement = DofTable::kPointsPerElement;
  static constexpr int knumNodesPerFace = DofTable::kNumNodesPerFace;
  static constexpr auto kFaceToElemDof = DofTable::kFaceToElemDof;
  static constexpr auto kFaceToElemDofAtDepth = DofTable::kFaceToElemDofAtDepth;
};

/// Acoustic DGsolver.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using DGsolverAcoustic =
    DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kAcoustic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_
