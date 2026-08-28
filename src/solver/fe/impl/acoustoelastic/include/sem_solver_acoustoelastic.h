#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_H_

#include <array>
#include <cmath>

#include "data_type.h"
#include "rhs_acoustoelastic.h"
#include "rhs_elastic.h"
#include "sem_enums.h"
#include "sem_solver.h"
#include "sem_solver_data.h"
#include "solver.h"
#include "wavefield_acoustoelastic.h"

namespace solver {
namespace fe {

/// Element belongs to the acoustic (fluid) domain.
static constexpr int kElementTypeAcoustic = 1;
/// Element belongs to the elastic (solid) domain.
static constexpr int kElementTypeElastic = 2;

/**
 * @brief Data structure passed to SEMsolverAcoustoElastic at each time step.
 *
 * Combines the acoustic wavefield, the elastic wavefield, and the acoustic
 * source term (source in the fluid domain only for this V1).
 */
struct SEMsolverDataAcoustoElastic : public Solver::DataStruct {
  /**
   * @param wavefield Combined acousto-elastic wavefield (p + ux/uy/uz).
   * @param rhs       Acoustic source term applied to the fluid domain.
   */
  SEMsolverDataAcoustoElastic(const WavefieldAcoustoElastic& wavefield, const RhsAcoustoElastic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const { return m_wavefield.getCurrentField(i); }

  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const { return m_wavefield.getPreviousField(i); }

  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const { return m_wavefield.getPrevPrevField(i); }

  void print() const override {
    m_wavefield.print();
    m_rhs.print();
  }

  /// Swap previous/current wavefields (call once per time step after
  /// computeOneStep).
  void swapWavefields() { m_wavefield.swap(); }

  WavefieldAcoustoElastic m_wavefield;  ///< Combined wavefield (p + u)
  RhsAcoustoElastic m_rhs;              ///< Acoustic source
};

/**
 * @brief Acousto-elastic coupled SEM solver (pressure–displacement
 * formulation, Komatitsch et al. 2000).
 *
 * Staggered explicit scheme: elastic step → A→E coupling → acoustic step →
 * E→A coupling.  Each sub-solver processes only its own elements via domain
 * masking.
 *
 * @tparam ORDER             Polynomial order of spectral elements.
 * @tparam INTEGRAL_TYPE     Quadrature/basis function type (Makutu kernels).
 * @tparam MESH_TYPE         Mesh implementation (ModelStruct or ModelUnstruct).
 * @tparam IS_MODEL_ON_NODES If true, material properties are stored on nodes.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
class SEMsolverAcoustoElastic : public Solver {
 public:
  using AcousticSolverType =
      SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kAcoustic>;
  using ElasticSolverType =
      SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kElastic>;
  using DataType = SEMsolverDataAcoustoElastic;

  SEMsolverAcoustoElastic() = default;
  ~SEMsolverAcoustoElastic() = default;

  // --- Solver interface ---

  int getNumComponents() const override { return 4; }  // p, ux, uy, uz

  /// @brief Returns the acoustic sub-solver mass matrix for DD synchronization.
  vectorReal& getMassMatrixAcoustic() override { return m_acoustic_solver_.getMassMatrixAcoustic(); }

  /// @brief Returns the elastic sub-solver mass matrix for DD synchronization.
  vectorReal& getMassMatrixElastic() override { return m_elastic_solver_.getMassMatrixElastic(); }

  vectorReal& getDampingMatrix(int c) override {
    if (c == 0) return m_acoustic_solver_.getDampingMatrix(0);
    return m_elastic_solver_.getDampingMatrix(c - 1);
  }

  vectorReal& getForceVector(int c) override {
    if (c == 0) return m_acoustic_solver_.getForceVector(0);
    return m_elastic_solver_.getForceVector(c - 1);
  }

  /// @brief Interface coupling coefficient c = int_Gamma phi n dGamma, one
  /// component per direction.  Assembled from the local acoustic element faces
  /// only, so a distributed driver must sum it over rank boundaries.
  vectorReal& getInterfaceCouplingCoeff(int c) override {
    if (c == 0) return m_coupling_coeff_x_;
    if (c == 1) return m_coupling_coeff_y_;
    return m_coupling_coeff_z_;
  }

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  /**
   * @brief Perform one coupled time step (serial / non-distributed mode).
   *
   * Implements the staggered explicit scheme:
   *   elastic forces → elastic update → acoustic forces → acoustic update,
   * with interface coupling applied between the two sub-steps.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override;

  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override;

  void updateSolutionForward(const float& dt, DataStruct& data) override;

  void updateSolutionBackward(const float& dt, DataStruct& data) override;

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;
  void computeDampingMatrix() override;

  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override;
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override {};

  void setAnisotropyType(model::AnisotropyType type) override { m_elastic_solver_.setAnisotropyType(type); }

  void setInterfacePropertyConvention(utils::enums::interfacePropertyConvention convention) override {
    interface_property_convention_ = convention;
  }

  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    m_acoustic_solver_.setSLSAttenuation(reference_frequencies, anelasticity_coefficients);
    m_elastic_solver_.setSLSAttenuation(reference_frequencies, anelasticity_coefficients);
  }

  // --- Accessors for diagnostics ---

  /// Number of acoustic elements detected in the mesh.
  int getNumAcousticElements() const { return num_acoustic_elements_; }

  /// Number of elastic elements detected in the mesh.
  int getNumElasticElements() const { return num_elastic_elements_; }

  /// Number of interface nodes (adjacent to both domains).
  int getNumInterfaceNodes() const { return num_interface_nodes_; }

  // GPU kernels must reside in public methods (CUDA extended lambda
  // constraint).

  /// @brief Identify interface nodes (adjacent to both domains).
  void TagNodes();

  /**
   * @brief Compute per-node interface coupling coefficients (area-weighted
   * outward normal integrated over each interface face).
   */
  void ComputeInterfaceCouplingCoefficients();

  /// @brief Snapshot u^{n-1} at the interface nodes before the elastic update
  /// overwrites it; needed by the elastic-to-acoustic coupling.
  void SaveInterfaceUnm1(const DataType& data);

  /**
   * @brief Apply acoustic→elastic coupling post-Verlet.
   * @param dt   Time step.
   * @param data Coupled solver data.
   */
  void ApplyCouplingAcousticToElastic(float dt, const DataType& data);

  /**
   * @brief Apply elastic→acoustic coupling post-Verlet.
   * @param data Coupled solver data.
   */
  void ApplyCouplingElasticToAcoustic(float dt, const DataType& data);

  /**
   * @brief Enforce the fluid/solid interface conditions on the two predictors.
   * @param dt   Time step.
   * @param data Coupled solver data, with both sub-domains already advanced.
   */
  void ApplyInterfaceCoupling(float dt, const DataType& data);

 private:
  AcousticSolverType m_acoustic_solver_;  ///< Acoustic sub-solver
  ElasticSolverType m_elastic_solver_;    ///< Elastic sub-solver

  MESH_TYPE m_mesh_;  ///< Local copy of the mesh

  /// Per-element type tag (kElementTypeAcoustic or kElementTypeElastic).
  vectorInt m_element_type_;

  /// Index map from global node index to interface node index (-1 if not).
  vectorInt m_interface_node_index_;

  /// Number of fluid–solid interface nodes.
  int n_interface_nodes_ = 0;
  /// Compact list of global interface node indices (size n_interface_nodes_).
  vectorInt m_interface_node_indices_;

  /// Area-weighted outward normal (solid→fluid) per node — X/Y/Z components.
  vectorReal m_coupling_coeff_x_;
  vectorReal m_coupling_coeff_y_;
  vectorReal m_coupling_coeff_z_;

  /// Elastic displacement at time n-1, stored for interface nodes only
  /// (size = n_interface_nodes_).  Allocated at the end of TagNodes.
  vectorReal m_ux_nm1_iface_;
  vectorReal m_uy_nm1_iface_;
  vectorReal m_uz_nm1_iface_;

  /// @brief One adjacent elastic element per interface node (size
  /// n_interface_nodes_).  Used to recover solid properties at interface nodes
  /// when IS_MODEL_ON_NODES is true.
  vectorInt m_interface_adj_elastic_elem_;

  /// @brief Solid material properties at interface nodes (size
  /// n_interface_nodes_).  Valid only when IS_MODEL_ON_NODES is true.
  vectorReal m_vp_solid_iface_;
  vectorReal m_vs_solid_iface_;
  vectorReal m_rho_solid_iface_;

  /// @brief Fluid material properties at interface nodes (size
  /// n_interface_nodes_).  Valid only when IS_MODEL_ON_NODES is true.
  vectorReal m_vp_fluid_iface_;
  vectorReal m_rho_fluid_iface_;

  /// @brief Set by the caller before computeFEInit; drives how TagNodes fills
  /// the solid side of the interface nodes.
  utils::enums::interfacePropertyConvention interface_property_convention_{
      utils::enums::interfacePropertyConvention::kFluidOnInterfaceNodes};

  int num_acoustic_elements_{0};  ///< Count of acoustic elements
  int num_elastic_elements_{0};   ///< Count of elastic elements
  int num_interface_nodes_{0};    ///< Count of interface nodes

  /// @brief Compact list of acoustic element indices (size
  /// num_acoustic_elements_).
  vectorInt acoustic_elem_list_;
  /// @brief Compact list of elastic element indices (size
  /// num_elastic_elements_).
  vectorInt elastic_elem_list_;

  int num_acoustic_nodes_{0};  ///< Count of acoustic-domain nodes
  int num_elastic_nodes_{0};   ///< Count of elastic-domain nodes

  /// @brief Compact list of acoustic-domain node indices (pure acoustic +
  /// interface, size num_acoustic_nodes_).
  vectorInt acoustic_node_list_;
  /// @brief Compact list of elastic-domain node indices (pure elastic +
  /// interface, size num_elastic_nodes_).
  vectorInt elastic_node_list_;

  /// Shear-modulus threshold below which an element is classified as acoustic.
  static constexpr float kMuTolerance = 1.0e-6f;

  /// @brief Classify each element as acoustic or elastic (mu < kMuTolerance).
  void TagElements();
};

}  // namespace fe
}  // namespace solver

#include "sem_solver_acoustoelastic_impl.h"

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_H_
