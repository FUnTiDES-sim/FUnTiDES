#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_H_

#include <array>

#include "data_type.h"
#include "dg-sem_physics_traits_acoustic.h"
#include "dg-sem_solver_data.h"
#include "dg_solver.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "sem_enums.h"
#include "sem_solver.h"
#include "solver.h"

namespace solver {
namespace fe {

/// Tag value of an element that belongs to the DG domain.
static constexpr int kElementTypeDG = 1;
/// Tag value of an element that belongs to the SEM domain.
static constexpr int kElementTypeSEM = 2;

/**
 * @brief Solver coupling a DG domain and a SEM domain on one mesh.
 *
 * Each time step applies, in this order: SEM to DG coupling, DG to SEM coupling, DG step,
 * SEM step. Each sub-solver only processes its own elements, given as a compact element list.
 *
 * @tparam ORDER             Polynomial order of elements.
 * @tparam INTEGRAL_TYPE     Quadrature/basis function type.
 * @tparam MESH_TYPE         Mesh type (model::ModelStruct or model::ModelUnstruct).
 * @tparam IS_MODEL_ON_NODES If true, material properties are stored on nodes.
 * @tparam PHYSICS           Physical model type.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
class DGSEMsolver : public Solver {
 public:
  using dgSolver = DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>;
  using semSolver = SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>;
  using DataType = DGSEMsolverData;

  DGSEMsolver() = default;
  ~DGSEMsolver() = default;

  /// Number of wavefield components.
  static constexpr int kNumFields = DGSEMPhysicsTraits::WavefieldType::kNumFields;

  int getNumComponents() const override { return kNumFields; }

  /// No-op: there is no global FE array to initialize.
  void initFEarrays() override {}

  /// No-op: sponge values are not handled by this solver.
  void initSpongeValues() override {}

  /// No-op: there is no global vector.
  void resetGlobalVectors(int numNodes) override {}

  /// No-op: there is no global mass matrix.
  void computeGlobalMassMatrix() override {}

  /// No-op: there is no global damping matrix.
  void computeDampingMatrix() override {}

  /// No-op: use computeOneStep().
  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override {}

  /// @throws std::runtime_error always: there is no global mass matrix.
  vectorReal& getMassMatrixAcoustic() override {
    throw std::runtime_error("getMassMatrixAcoustic not implemented for DG");
  }

  /// @throws std::runtime_error always: there is no global mass matrix.
  vectorReal& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG-SEM coupling");
  }

  /// @throws std::runtime_error always: there is no global damping matrix.
  vectorReal& getDampingMatrix(int c) override {
    throw std::runtime_error("getDampingMatrix not implemented for DG");
  }

  /// @throws std::runtime_error always: there is no global force vector.
  vectorReal& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
  }

  /// No-op: use computeOneStep().
  void updateSolutionForward(const float& dt, DataStruct& data) override {}

  /// No-op: there is no backward propagation.
  void updateSolutionBackward(const float& dt, DataStruct& data) override {}

  /// No-op.
  void setAnisotropyType(model::AnisotropyType type) override {
    // TODO: anisotropy is not supported by the coupled solver yet.
  }

  /// Set the Z coordinate of the DG-SEM interface used by TagElements().
  void setZBoundary(float z) override { DG_SEM_interface_z_ = z; }

  /**
   * @brief Provide the DG/SEM element split directly, bypassing the Z-threshold split of
   * TagElements().
   *
   * Must be called before computeFEInit(). If not called, or if the size of @p tags differs
   * from the mesh element count, TagElements() uses the Z-threshold split. The threshold probes
   * the deformed Z coordinate of a single node, so it only cuts the intended plane on a flat
   * mesh; a caller that knows the split (for example from a layer index) should use this method.
   *
   * @param[in] tags Element type (kElementTypeDG or kElementTypeSEM), one entry per mesh element.
   */
  void setElementTags(const vectorInt& tags) override { m_external_element_type_ = tags; }

  /// No-op.
  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    // TODO: SLS attenuation is not supported by the coupled solver yet.
  }

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  void allocateFEarrays() override;

  /// @brief Identify interface nodes (adjacent to both domains).
  void TagNodes();

  /// @brief Classify each element as DG or SEM.
  void TagElements();

  /**
   * @brief Perform one coupled time step.
   *
   * Applies the interface coupling, then the DG step and the SEM step.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override;

  /**
   * @brief Compute the SIPG interface flux contribution on both sides (DG and SEM).
   *
   * Reads the current-step pressure from both domains. The contribution is accumulated into the
   * DG local stiffness buffer (m_stiff_local_) and into the SEM workVectorsGlobal_[0].
   *
   * @param[in] data Coupled solver data.
   */
  void ApplyCoupling(const DataType& data);

  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override;
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override;

  /// @return Number of DG elements detected in the mesh.
  int getNumDGElements() const { return num_DG_elements_; }

  /// @return Number of SEM elements detected in the mesh.
  int getNumSEmElements() const { return num_SEm_elements_; }

  /// @return Number of interface faces (adjacent to both domains).
  int getNumInterfaceFaces() const { return num_interface_faces_; }

 private:
  dgSolver m_DG_solver_;
  semSolver m_SEm_solver_;

  MESH_TYPE m_mesh_;
  /// Shared with m_DG_solver_, see DGsolver::setFaceConnectivity.
  model::FaceConnectivityUnstruct<float, int, ORDER> m_face_connectivity_;
  static constexpr int knumNodesPerFace = (ORDER + 1) * (ORDER + 1);

  /// Element type (kElementTypeDG or kElementTypeSEM), one entry per mesh element.
  vectorInt m_element_type_;

  int num_interface_faces_{0};
  /// Global indices of the interface faces, size num_interface_faces_.
  vectorInt m_interface_face_indices_;

  int num_DG_elements_{0};
  int num_SEm_elements_{0};
  /// Indices of the DG elements, size num_DG_elements_.
  vectorInt DG_elem_list_;
  /// Indices of the SEM elements, size num_SEm_elements_.
  vectorInt SEm_elem_list_;

  int num_SEm_nodes_{0};
  /// Indices of the SEM-domain nodes (pure SEM plus interface), size num_SEm_nodes_.
  vectorInt SEm_node_list_;

  int m_n_DG_interior_faces_{0};  ///< Number of DG-DG faces, interface faces excluded.
  /// Indices of the DG-DG faces, interface faces excluded, size m_n_DG_interior_faces_.
  vectorInt m_DG_interior_face_list_;

  /// @brief Build m_DG_interior_face_list_ from all DG faces minus the interface faces.
  void BuildDGInteriorFaceList();

  float DG_SEM_interface_z_ = 1000.f;  ///< Z coordinate of the DG-SEM interface.
  ///< Element types given by the caller through setElementTags(). Empty unless set.
  vectorInt m_external_element_type_;
  ///< SIPG penalty factor; must stay equal to the value held by the DG sub-solver.
  real_t m_penalty_factor_ = 12.0f;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_H_
