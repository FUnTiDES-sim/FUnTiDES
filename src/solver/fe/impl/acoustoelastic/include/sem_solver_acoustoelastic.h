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

static constexpr int kElementTypeAcoustic = 1;  ///< Element tag: acoustic (fluid) domain.
static constexpr int kElementTypeElastic = 2;   ///< Element tag: elastic (solid) domain.

/**
 * @brief Data passed to SEMsolverAcoustoElastic at each time step.
 *
 * Bundles the combined acousto-elastic wavefield and the source term. The source
 * is applied to the fluid domain only.
 */
struct SEMsolverDataAcoustoElastic : public Solver::DataStruct {
  /**
   * @param wavefield Combined wavefield (pressure and displacement components).
   * @param rhs       Source term.
   */
  SEMsolverDataAcoustoElastic(const WavefieldAcoustoElastic& wavefield, const RhsAcoustoElastic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  /// @brief Field component i at the current time step.
  PROXY_HOST_DEVICE
  vectorReal getCurrentField(int i) const { return m_wavefield.getCurrentField(i); }

  /// @brief Field component i at the previous time step.
  PROXY_HOST_DEVICE
  vectorReal getPreviousField(int i) const { return m_wavefield.getPreviousField(i); }

  /// @brief Field component i two time steps back.
  PROXY_HOST_DEVICE
  vectorReal getPrevPrevField(int i) const { return m_wavefield.getPrevPrevField(i); }

  void print() const override {
    m_wavefield.print();
    m_rhs.print();
  }

  /// @brief Rotate the wavefield time levels; call once per time step after computeOneStep.
  void swapWavefields() { m_wavefield.swap(); }

  WavefieldAcoustoElastic m_wavefield;  ///< Combined wavefield.
  RhsAcoustoElastic m_rhs;              ///< Source term.
};

/**
 * @brief Coupled acoustic/elastic SEM solver (pressure-displacement formulation,
 * Komatitsch et al. 2000).
 *
 * One time step is a staggered explicit scheme: elastic step, acoustic-to-elastic
 * coupling, acoustic step, elastic-to-acoustic coupling. Each sub-solver only
 * processes the elements of its own domain.
 *
 * @tparam ORDER             Polynomial order of the spectral elements.
 * @tparam INTEGRAL_TYPE     Quadrature and basis function type.
 * @tparam MESH_TYPE         Mesh type.
 * @tparam IS_MODEL_ON_NODES True if material properties are stored on nodes, false if on elements.
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

  /// @brief Number of field components: pressure, then ux, uy, uz.
  int getNumComponents() const override { return 4; }

  /// @brief Mass matrix of the acoustic domain, exposed for domain-decomposition synchronization.
  vectorReal& getMassMatrixAcoustic() override { return m_acoustic_solver_.getMassMatrixAcoustic(); }

  /// @brief Mass matrix of the elastic domain, exposed for domain-decomposition synchronization.
  vectorReal& getMassMatrixElastic() override { return m_elastic_solver_.getMassMatrixElastic(); }

  /// @brief Damping matrix of component c (0 is pressure, 1 to 3 are the displacement components).
  vectorReal& getDampingMatrix(int c) override {
    if (c == 0) return m_acoustic_solver_.getDampingMatrix(0);
    return m_elastic_solver_.getDampingMatrix(c - 1);
  }

  /// @brief Force vector of component c (0 is pressure, 1 to 3 are the displacement components).
  vectorReal& getForceVector(int c) override {
    if (c == 0) return m_acoustic_solver_.getForceVector(0);
    return m_elastic_solver_.getForceVector(c - 1);
  }

  /**
   * @brief Interface coupling coefficient c = int_Gamma phi n dGamma, direction c (0=x, 1=y, 2=z).
   *
   * Assembled from the local acoustic element faces only, so a distributed
   * driver must sum it over rank boundaries.
   */
  vectorReal& getInterfaceCouplingCoeff(int c) override {
    if (c == 0) return m_coupling_coeff_x_;
    if (c == 1) return m_coupling_coeff_y_;
    return m_coupling_coeff_z_;
  }

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  /**
   * @brief Perform one coupled time step in serial (non-distributed) mode.
   *
   * Elastic forces, elastic update, acoustic forces, acoustic update, with the
   * interface coupling applied between the two updates.
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

  /// @brief Select how the solid side of the interface nodes is filled; call before computeFEInit.
  void setInterfacePropertyConvention(utils::enums::interfacePropertyConvention convention) override {
    interface_property_convention_ = convention;
  }

  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    m_acoustic_solver_.setSLSAttenuation(reference_frequencies, anelasticity_coefficients);
    m_elastic_solver_.setSLSAttenuation(reference_frequencies, anelasticity_coefficients);
  }

  /// @brief Number of acoustic elements detected in the mesh.
  int getNumAcousticElements() const { return num_acoustic_elements_; }

  /// @brief Number of elastic elements detected in the mesh.
  int getNumElasticElements() const { return num_elastic_elements_; }

  /// @brief Number of interface nodes (nodes adjacent to both domains).
  int getNumInterfaceNodes() const { return num_interface_nodes_; }

  // The methods below are public because GPU kernels (extended lambdas) cannot
  // live in private or protected methods.

  /// @brief Identify the interface nodes and the node lists of each domain.
  void TagNodes();

  /// @brief Compute the per-node interface coupling coefficients (area-weighted
  /// outward normal integrated over each interface face).
  void ComputeInterfaceCouplingCoefficients();

  /// @brief Save the elastic displacement at time n-1 on the interface nodes,
  /// before the elastic update overwrites it; needed by the elastic-to-acoustic coupling.
  void SaveInterfaceUnm1(const DataType& data);

  /**
   * @brief Apply the acoustic-to-elastic coupling after the Verlet update.
   * @param dt   Time step.
   * @param data Coupled solver data.
   */
  void ApplyCouplingAcousticToElastic(float dt, const DataType& data);

  /**
   * @brief Apply the elastic-to-acoustic coupling after the Verlet update.
   * @param dt   Time step.
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
  AcousticSolverType m_acoustic_solver_;  ///< Solver of the acoustic domain.
  ElasticSolverType m_elastic_solver_;    ///< Solver of the elastic domain.

  MESH_TYPE m_mesh_;  ///< Copy of the mesh.

  /// Per-element domain tag (kElementTypeAcoustic or kElementTypeElastic).
  vectorInt m_element_type_;

  /// Global node index to interface node index, -1 if the node is not on the interface.
  vectorInt m_interface_node_index_;

  int n_interface_nodes_ = 0;  ///< Number of fluid/solid interface nodes.
  /// Global indices of the interface nodes, size n_interface_nodes_.
  vectorInt m_interface_node_indices_;

  /// Area-weighted outward normal (solid to fluid), x component.
  vectorReal m_coupling_coeff_x_;
  /// Area-weighted outward normal (solid to fluid), y component.
  vectorReal m_coupling_coeff_y_;
  /// Area-weighted outward normal (solid to fluid), z component.
  vectorReal m_coupling_coeff_z_;

  /// Elastic displacement at time n-1 on the interface nodes, size n_interface_nodes_;
  /// allocated at the end of TagNodes. One vector per direction (x, y, z).
  vectorReal m_ux_nm1_iface_;
  vectorReal m_uy_nm1_iface_;
  vectorReal m_uz_nm1_iface_;

  /// One adjacent elastic element per interface node, size n_interface_nodes_. Used to
  /// recover the solid properties at interface nodes when IS_MODEL_ON_NODES is true.
  vectorInt m_interface_adj_elastic_elem_;

  /// Solid vp, vs and rho at the interface nodes, size n_interface_nodes_.
  /// Valid only when IS_MODEL_ON_NODES is true.
  vectorReal m_vp_solid_iface_;
  vectorReal m_vs_solid_iface_;
  vectorReal m_rho_solid_iface_;

  /// Fluid vp and rho at the interface nodes, size n_interface_nodes_.
  /// Valid only when IS_MODEL_ON_NODES is true.
  vectorReal m_vp_fluid_iface_;
  vectorReal m_rho_fluid_iface_;

  /// Drives how TagNodes fills the solid side of the interface nodes; must be set before computeFEInit.
  utils::enums::interfacePropertyConvention interface_property_convention_{
      utils::enums::interfacePropertyConvention::kFluidOnInterfaceNodes};

  int num_acoustic_elements_{0};  ///< Number of acoustic elements.
  int num_elastic_elements_{0};   ///< Number of elastic elements.
  int num_interface_nodes_{0};    ///< Number of interface nodes.

  /// Indices of the acoustic elements, size num_acoustic_elements_.
  vectorInt acoustic_elem_list_;
  /// Indices of the elastic elements, size num_elastic_elements_.
  vectorInt elastic_elem_list_;

  int num_acoustic_nodes_{0};  ///< Number of acoustic-domain nodes.
  int num_elastic_nodes_{0};   ///< Number of elastic-domain nodes.

  /// Acoustic-domain node indices (purely acoustic and interface), size num_acoustic_nodes_.
  vectorInt acoustic_node_list_;
  /// Elastic-domain node indices (purely elastic and interface), size num_elastic_nodes_.
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
