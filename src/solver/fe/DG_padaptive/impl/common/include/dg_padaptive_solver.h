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

static constexpr int kElementTypePMin = 1;  ///< Element tag: element belongs to the pMin order domain.
static constexpr int kElementTypePMax = 2;  ///< Element tag: element belongs to the pMax order domain.

/**
 * @brief Acoustic DG solver with two polynomial orders (pMin and pMax) on one mesh.
 *
 * Each time step couples the two domains through a SIPG interface flux, then advances the pMin
 * and pMax sub-solvers in turn. Each sub-solver only processes its own list of elements.
 *
 * @tparam ORDER_MIN          Polynomial order pMin of the low-order elements.
 * @tparam ORDER_MAX          Polynomial order pMax of the high-order elements.
 * @tparam INTEGRAL_SELECTOR  Template mapping (order, IMPL_TAG) to the quadrature/basis type.
 * @tparam IMPL_TAG           Implementation tag passed to INTEGRAL_SELECTOR.
 * @tparam MESH_TYPE          Mesh implementation.
 * @tparam IS_MODEL_ON_NODES  If true, material properties are stored on nodes.
 * @tparam PHYSICS            Physical model type.
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

  static constexpr int kNumFields = DGPAdaptivePhysicsTraits::WavefieldType::kNumFields;  ///< Fields per wavefield.

  /// Number of 1D dofs per direction in the pMin domain. Class constant because nvcc rejects a
  /// kernel-local constexpr used as a host identifier inside a device lambda.
  static constexpr int kNumDofs1dMin = ORDER_MIN + 1;
  /// Number of 1D dofs per direction in the pMax domain. Same rationale as kNumDofs1dMin.
  static constexpr int kNumDofs1dMax = ORDER_MAX + 1;

  /// @return Number of wavefield components.
  int getNumComponents() const override { return kNumFields; }

  /// No-op: this solver has no global finite-element arrays.
  void initFEarrays() override {}

  /// No-op: sponge values are not used by this solver.
  void initSpongeValues() override {}

  /// No-op: this solver has no global vectors.
  void resetGlobalVectors(int numNodes) override {}

  /// No-op: this solver has no global mass matrix.
  void computeGlobalMassMatrix() override {}

  /// No-op: this solver has no global damping matrix.
  void computeDampingMatrix() override {}

  /// No-op: forces are applied inside computeOneStep().
  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override {}

  /// @throws std::runtime_error Always: there is no global mass matrix.
  vectorReal& getMassMatrixAcoustic() override {
    throw std::runtime_error("getMassMatrixAcoustic not implemented for DG");
  }

  /// @throws std::runtime_error Always: there is no global mass matrix.
  vectorReal& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG");
  }

  /// @throws std::runtime_error Always: there is no global damping matrix.
  vectorReal& getDampingMatrix(int c) override { throw std::runtime_error("getDampingMatrix not implemented for DG"); }

  /// @throws std::runtime_error Always: there is no global force vector.
  vectorReal& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
  }

  /// No-op: the update is done by computeOneStep().
  void updateSolutionForward(const float& dt, DataStruct& data) override {}

  /// No-op: backward propagation is not supported.
  void updateSolutionBackward(const float& dt, DataStruct& data) override {}

  /// No-op: anisotropy is not supported.
  void setAnisotropyType(model::AnisotropyType type) override {
    // TODO: anisotropy is not supported by the DG solvers.
  }

  /// @param z Z coordinate of the pMin/pMax interface, used by TagElements() when no element tags were given.
  void setZBoundary(float z) override { pAdaptive_interface_z_ = z; }

  /**
   * @brief Provide the pMin/pMax element split directly instead of the Z-threshold split.
   *
   * Must be called before computeFEInit(), which calls TagElements(). If the size of the tags differs
   * from the mesh element count (or if this is never called), TagElements() falls back to the
   * Z-threshold split. The Z-threshold probes the deformed coordinate of a single node per
   * element, so it only cuts the intended plane on a flat mesh.
   *
   * @param[in] tags Per-element type (kElementTypePMin or kElementTypePMax), one entry per mesh element.
   */
  void setElementTags(const vectorInt& tags) override { m_external_element_type_ = tags; }

  /// No-op: attenuation is not supported.
  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    // TODO: SLS attenuation is not supported by the DG solvers.
  }

  /**
   * @brief Tag elements and interface nodes, build the sub-solvers and the interface data.
   *
   * @param[in,out] mesh Mesh, must be of type MESH_TYPE.
   * @todo VERIFY: are sponge_size, surface_sponge and taper_delta ignored on purpose?
   */
  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  /// Allocate the arrays of the sub-solvers and of the interface coupling.
  void allocateFEarrays() override;

  /// Identify the interface nodes (adjacent to both domains).
  void TagNodes();

  /// Classify each element as pMin or pMax order and build the element lists.
  void TagElements();

  /// Fill the 1D order-raising matrix m_p1d_projection_, from which the mortar projection is built.
  void ComputeMortarProjection();

  /**
   * @brief Raise the pMin pressure to ORDER_MAX on every interface-adjacent pMin element, and
   * zero the matching interface stiffness accumulator.
   *
   * The raise is exact: P_ORDER_MIN is a subspace of P_ORDER_MAX. The whole element is raised,
   * not only the face layer, because the SIPG consistency term reads the normal derivative,
   * which lives on the depth line behind each face dof.
   *
   * @param[in] data Coupled solver data.
   */
  void ProlongPMinField(const DataType& data);

  /// Restrict the interface stiffness from the fictitious ORDER_MAX grid onto the real pMin dofs
  /// and accumulate it into the pMin sub-solver. Adjoint of ProlongPMinField().
  void RestrictPMinStiff();

  /**
   * @brief Perform one coupled time step (non-distributed mode).
   *
   * Applies the interface coupling, then advances the pMin domain, then the pMax domain.
   *
   * @param[in] dt         Time step.
   * @param[in] timeSample Index of the current time sample.
   * @param[in,out] data   Must be a DataType.
   * @throws std::bad_cast If data is not a DataType.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override;

  /**
   * @brief Compute the SIPG interface flux between the pMax and pMin domains.
   *
   * Reads the pressure of the current step from both domains and accumulates the result into the
   * local stiffness of the DG sub-solvers, which applyVerlet consumes.
   *
   * @param[in] data Coupled solver data.
   */
  void ApplyCoupling(const DataType& data);

  /// No-op for a vectorReal field.
  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override {};
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override;

  /// @return Number of pMin elements in the mesh.
  int getNumPMinElements() const { return num_pMin_elements_; }

  /// @return Number of pMax elements in the mesh.
  int getNumPMaxElements() const { return num_pMax_elements_; }

  /// @return Number of interface faces (adjacent to both domains).
  int getNumInterfaceFaces() const { return num_interface_faces_; }

 private:
  pMinSolver m_pMin_solver_;  ///< Sub-solver of the pMin domain.
  pMaxSolver m_pMax_solver_;  ///< Sub-solver of the pMax domain.

  vectorInt order_list;

  MESH_TYPE m_mesh_;  ///< Local copy of the mesh, built at the highest order.
  model::FaceConnectivityUnstruct<float, int, ORDER_MAX> m_face_connectivity_;  ///< Face connectivity at ORDER_MAX.

  /// 1D order-raising matrix, m_p1d_projection_(k, m) = phi^pMin_m(xi^pMax_k). The mortar
  /// projection of ApplyCoupling is its threefold tensor product.
  arrayReal m_p1d_projection_;

  /// Compact list of the pMin elements touching the interface, the only ones the coupling reads
  /// at ORDER_MAX resolution. Sized by the interface surface, not the volume.
  vectorInt m_iface_pMin_elem_list_;
  /// Row of the prolonged arrays owned by each pMin element, -1 away from the interface.
  vectorInt m_pMin_elem_to_slot_;
  int m_n_iface_pMin_elements_{0};  ///< Number of interface-adjacent pMin elements.

  /// pMin pressure raised to ORDER_MAX, one row per interface-adjacent pMin element.
  arrayReal m_pMin_prolonged_field_;
  /// Interface stiffness accumulated on the fictitious ORDER_MAX grid, restricted onto the real
  /// pMin dofs at the end of every step.
  arrayReal m_pMin_prolonged_stiff_;

  vectorInt m_element_type_;  ///< Per-element type tag (kElementTypePMin or kElementTypePMax).

  int num_interface_faces_{0};  ///< Number of interface faces.
  /// Global indices of the interface faces, size num_interface_faces_.
  vectorInt m_interface_face_indices_;

  int num_pMin_elements_{0};  ///< Number of pMin elements.
  int num_pMax_elements_{0};  ///< Number of pMax elements.
  vectorInt pMin_elem_list_;  ///< Indices of the pMin elements, size num_pMin_elements_.
  vectorInt pMax_elem_list_;  ///< Indices of the pMax elements, size num_pMax_elements_.

  int m_n_pMin_interior_faces_{0};  ///< Number of pMin-pMin interior faces (excludes interface faces).
  int m_n_pMax_interior_faces_{0};  ///< Number of pMax-pMax interior faces (excludes interface faces).
  /// Indices of the pMin-pMin interior faces.
  vectorInt m_pMin_interior_face_list_;
  /// Indices of the pMax-pMax interior faces.
  vectorInt m_pMax_interior_face_list_;

  /// Build m_pMin_interior_face_list_ and m_pMax_interior_face_list_ from the element types:
  /// all faces minus the interface faces.
  void BuildInteriorFaceLists();

  /// Build m_iface_pMin_elem_list_ and m_pMin_elem_to_slot_ and allocate the prolonged arrays.
  /// Requires the interface face list.
  void BuildInterfaceElementList();

  float pAdaptive_interface_z_ = 1000.f;  ///< Z coordinate of the pMin/pMax interface (Z-threshold split).
  /// Caller-provided per-element type tags, see setElementTags(). Empty unless set.
  vectorInt m_external_element_type_;
  /// SIPG penalty factor of the interface coupling. Overwritten in computeFEInit() with the
  /// penalty factor of the sub-solvers.
  real_t m_penalty_factor_ = 12.0f;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_H_
