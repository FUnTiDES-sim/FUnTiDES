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

/// Element belongs to the DG domain.
static constexpr int kElementTypeDG = 1;
/// Element belongs to the SEM domain.
static constexpr int kElementTypeSEM = 2;

/**
 * @brief DG-SEm coupled solver.
 *
 * Staggered explicit scheme: SEM→DG coupling → DG→SEM coupling → DG step → SEM step.
 * Each sub-solver processes only its own elements via a list of elements.
 *
 * @tparam ORDER             Polynomial order of elements.
 * @tparam INTEGRAL_TYPE     Quadrature/basis function type (Makutu kernels).
 * @tparam MESH_TYPE         Mesh implementation (ModelStruct or ModelUnstruct).
 * @tparam IS_MODEL_ON_NODES If true, material properties are stored on nodes.
 * @tparam PHYSICS           Physical model type (Acoustic).
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

  // --- Solver interface ---

  static constexpr int kNumFields = DGSEMPhysicsTraits::WavefieldType::kNumFields;

  int getNumComponents() const override { return kNumFields; }

  // --- Mandatory overrides for Solver interface ---

  void initFEarrays() override {
    // Here for retrocompatibility
  }

  void initSpongeValues() override {
    // Here for retrocompatibility
  }

  void resetGlobalVectors(int numNodes) override {
    // Here for retrocompatibility, no global vector for DG-SEM
  }

  void computeGlobalMassMatrix() override {
    // Here for retrocompatibility, no global mass matrix for DG-SEM
  }

  void computeDampingMatrix() override {
    // Here for retrocompatibility, no global damping matrix for DG-SEM
  }

  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override {
    // Here for retrocompatibility
  }

  vectorReal& getMassMatrixAcoustic() override {
    throw std::runtime_error("getMassMatrixAcoustic not implemented for DG");
    // Maybe return SEM mass matrix
  }

  vectorReal& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG-SEM coupling");
  }

  vectorReal& getDampingMatrix(int c) override {
    throw std::runtime_error("getDampingMatrix not implemented for DG");
    // Maybe return SEM damping matrix
  }

  vectorReal& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
    // Maybe return SEM force vector
  }

  void updateSolutionForward(const float& dt, DataStruct& data) override {
    // Here for retrocompatibility
  }

  void updateSolutionBackward(const float& dt, DataStruct& data) override {
    // Here for retrocompatibility
  }

  void setAnisotropyType(model::AnisotropyType type) override {
    // TODO: Implement anisotropy setting
  }

  void setZBoundary(float z) override { DG_SEM_interface_z_ = z; }

  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    // TODO: Implement SLS attenuation setting
  }

  // --- Core solver methods ---

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  void allocateFEarrays() override;

  /// @brief Identify interface nodes (adjacent to both domains).
  void TagNodes();

  /// @brief Classify each element as SEm or DG.
  void TagElements();

  /**
   * @brief Perform one coupled time step (serial / non-distributed mode).
   *
   * Implements the staggered explicit scheme
   * with interface coupling applied between the two sub-steps.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override;

  /**
   * @brief Compute SIPG interface flux contribution on both side (DG and SEM).
   *
   * Reads p^n from both domains (no temporal lag). Accumulates into DG m_stiff_local_ and SEM workVectorsGlobal_[0],
   * consumed by applyVerlet.
   * @param data Coupled solver data.
   */
  void ApplyCoupling(const DataType& data);

  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override;
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override;

  // --- Accessors for diagnostics ---

  /// Number of DG elements detected in the mesh.
  int getNumDGElements() const { return num_DG_elements_; }

  /// Number of SEM elements detected in the mesh.
  int getNumSEmElements() const { return num_SEm_elements_; }

  /// Number of interface faces (adjacent to both domains).
  int getNumInterfaceFaces() const { return num_interface_faces_; }

 private:
  dgSolver m_DG_solver_;
  semSolver m_SEm_solver_;

  MESH_TYPE m_mesh_;
  /// Shared with m_DG_solver_ (see setFaceConnectivity): face id lists built
  /// here are consumed by the DG solver's face kernels, so both must index the
  /// same numbering. ORDER must match m_DG_solver_'s for the types to line up.
  model::FaceConnectivityUnstruct<float, int, ORDER> m_face_connectivity_;
  static constexpr int knumNodesPerFace = (ORDER + 1) * (ORDER + 1);

  /// Per-element type tag (kElementTypeDG or kElementTypeSEM).
  vectorInt m_element_type_;

  int num_interface_faces_{0};
  /// Compact list of global interface face indices (size num_interface_faces_).
  vectorInt m_interface_face_indices_;

  int num_DG_elements_{0};
  int num_SEm_elements_{0};
  /// Compact list of DG element indices (size num_DG_elements_).
  vectorInt DG_elem_list_;
  /// Compact list of SEm element indices (size num_SEm_elements_).
  vectorInt SEm_elem_list_;

  int num_SEm_nodes_{0};
  /// Compact list of SEm-domain node indices (pure SEm + interface, size num_SEm_nodes_).
  vectorInt SEm_node_list_;

  int m_n_DG_interior_faces_{0};  ///< Excludes DG-SEM interface faces
  /// Compact list of DG interior face indices (DG-DG faces only, excludes interface).
  vectorInt m_DG_interior_face_list_;

  /// @brief Build m_DG_interior_face_list_ from all DG faces minus interface faces.
  void BuildDGInteriorFaceList();

  float DG_SEM_interface_z_ = 1000.f;  ///< Z coordinate of the DG-SEM interface
  real_t m_penalty_factor_ = 12.0f;    ///< SIPG penalty; kept in sync with the DG sub-solver's own
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_H_
