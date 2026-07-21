#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_H_

#include <array>
#include <cmath>

#include "dg_padaptive_physics_traits_acoustic.h"
#include "dg_padaptive_solver_data.h"
#include "dg_solver.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "sem_enums.h"
#include "solver.h"

namespace solver {
namespace fe {

/// Element belongs to the pMin order domain.
static constexpr int kElementTypePMin = 1;
/// Element belongs to the pMax order domain.
static constexpr int kElementTypePMax = 2;

/**
 * @brief p-adaptive DG solver.
 *
 * Staggered explicit scheme: pMax→pMin coupling and pMin→pMax coupling → pMin step → pMax step.
 * Each sub-solver processes only its own elements via a list of elements.
 *
 * @tparam ORDER_MIN          Polynomial order pMin of elements.
 * @tparam ORDER_MAX          Polynomial order pMax of elements.
 * @tparam INTEGRAL_TYPE      Quadrature/basis function type (Makutu kernels).
 * @tparam MESH_TYPE          Mesh implementation (ModelStruct or ModelUnstruct).
 * @tparam IS_MODEL_ON_NODES  If true, material properties are stored on nodes.
 * @tparam PHYSICS            Physical model type (Acoustic).
 */
template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
class DGPAdaptiveSolver : public Solver {
 public:
  using INTEGRAL_TYPE_MIN = typename INTEGRAL_SELECTOR<ORDER_MIN, IMPL_TAG>::type;
  using INTEGRAL_TYPE_MAX = typename INTEGRAL_SELECTOR<ORDER_MAX, IMPL_TAG>::type;

  using pMinSolver = DGsolver<ORDER_MIN, INTEGRAL_TYPE_MIN, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>;
  using pMaxSolver = DGsolver<ORDER_MAX, INTEGRAL_TYPE_MAX, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>;
  using DataType = DGPAdaptiveSolverData;

  DGPAdaptiveSolver() = default;
  ~DGPAdaptiveSolver() = default;

  // --- Solver interface ---

  static constexpr int kNumFields = DGPAdaptivePhysicsTraits::WavefieldType::kNumFields;

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
  }

  vectorReal& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG");
  }

  vectorReal& getDampingMatrix(int c) override { throw std::runtime_error("getDampingMatrix not implemented for DG"); }

  vectorReal& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
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

  void setZBoundary(float z) override { pAdaptive_interface_z_ = z; }

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

  /// @brief Classify each element as pMin or pMax order.
  void TagElements();

  /// @brief Fill the m_mortar_projection matrix.
  void ComputeMortarProjection();

  /**
   * @brief Perform one coupled time step (serial / non-distributed mode).
   *
   * Implements the staggered explicit scheme
   * with interface coupling applied between the two sub-steps.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override;

  /**
   * @brief Compute SIPG interface flux contribution between pMax and pMin.
   *
   * Reads p^n from both domains (no temporal lag). Accumulates into DG m_stiff_local_,
   * consumed by applyVerlet.
   * @param data Coupled solver data.
   */
  void ApplyCoupling(const DataType& data);

  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override {};
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override;

  // --- Accessors for diagnostics ---

  /// Number of pMin elements detected in the mesh.
  int getNumPMinElements() const { return num_pMin_elements_; }

  /// Number of pMax elements detected in the mesh.
  int getNumPMaxElements() const { return num_pMax_elements_; }

  /// Number of interface faces (adjacent to both domains).
  int getNumInterfaceFaces() const { return num_interface_faces_; }

 private:
  pMinSolver m_pMin_solver_;  ///< pMin sub-solver
  pMaxSolver m_pMax_solver_;  ///< pMax sub-solver

  vectorInt order_list;

  MESH_TYPE m_mesh_;  ///< Local copy of the mesh built on the highest order
  model::FaceConnectivityUnstruct<float, int> m_face_connectivity_;  // pMax face connectivity structure

  /// pMin projection matrix: pressure field → mortar space
  arrayReal m_mortar_projection;

  /// Per-element type tag (kElementTypePMin or kElementTypePMax).
  vectorInt m_element_type_;

  int num_interface_faces_{0};  ///< Count of interface faces
  /// Compact list of global interface face indices (size n_interface_faces_).
  vectorInt m_interface_face_indices_;

  int num_pMin_elements_{0};  ///< Count of pMin elements
  int num_pMax_elements_{0};  ///< Count of pMax elements
  /// @brief Compact list of pMin element indices (size
  /// num_pMin_elements_).
  vectorInt pMin_elem_list_;
  /// @brief Compact list of pMax element indices (size
  /// num_pMax_elements_).
  vectorInt pMax_elem_list_;

  int m_n_pMin_interior_faces_{0};  ///< Count of pMin interior faces (excludes pMin-pMax interface faces)
  int m_n_pMax_interior_faces_{0};  ///< Count of pMax interior faces (excludes pMin-pMax interface faces)
  /// @brief Compact list of pMin interior face indices (pMin-pMin faces only, excludes interface).
  vectorInt m_pMin_interior_face_list_;
  /// @brief Compact list of pMax interior face indices (pMax-pMax faces only, excludes interface).
  vectorInt m_pMax_interior_face_list_;

  /// @brief Build based on kElementType, m_pMin_interior_face_list_
  /// and m_pMax_interior_face_list_ from all faces minus interface faces.
  void BuildInteriorFaceLists();

  float pAdaptive_interface_z_ = 1000.f;  ///< Z coordinate of the pMin-pMax interface
  /// @brief SIPG penalty factor for interface coupling (matches DG internal penalty of both solvers).
  real_t m_penalty_factor_ = 75.0f;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_H_
