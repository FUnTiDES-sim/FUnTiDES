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

namespace solver
{
namespace fe
{

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
struct SEMsolverDataAcoustoElastic : public Solver::DataStruct
{
  /**
   * @param wavefield Combined acousto-elastic wavefield (p + ux/uy/uz).
   * @param rhs       Acoustic source term applied to the fluid domain.
   */
  SEMsolverDataAcoustoElastic(const WavefieldAcoustoElastic& wavefield,
                              const RhsAcoustoElastic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs)
  {
  }

  void print() const override
  {
    m_wavefield.print();
    m_rhs.print();
  }

  /// Swap previous/current wavefields (call once per time step after computeOneStep).
  void swapWavefields() { m_wavefield.swap(); }

  WavefieldAcoustoElastic m_wavefield;  ///< Combined wavefield (p + u)
  RhsAcoustoElastic m_rhs;             ///< Acoustic source
};

/**
 * @brief Acousto-elastic coupled SEM solver (pressure–displacement formulation).
 *
 * Implements the elasto-acoustic coupling of Komatitsch et al. (2000) using
 * an explicit staggered time-stepping scheme:
 *
 *   1. Elastic forces step (includes acoustic pressure traction at interface).
 *   2. Elastic solution update.
 *   3. Acoustic forces step (includes elastic normal acceleration at interface).
 *   4. Acoustic solution update.
 *
 * The two sub-solvers (acoustic and elastic) each process only their own
 * elements via element masking. Interface nodes carry pre-computed outward
 * unit normals (solid→fluid direction).
 *
 * Mass matrix and damping matrix are re-computed with domain masking after the
 * sub-solver initialisation (sub-solvers accumulate over all elements; the
 * override restricts each sub-solver to its own domain).
 *
 * Solution update is delegated to the sub-solvers' updateFields(), which
 * correctly skips nodes whose mass is zero (i.e. nodes belonging to the
 * other domain).
 *
 * @tparam ORDER             Polynomial order of spectral elements (1, 2, or 3).
 * @tparam INTEGRAL_TYPE     Quadrature/basis function type (Makutu kernels).
 * @tparam MESH_TYPE         Mesh implementation (ModelStruct or ModelUnstruct).
 * @tparam IS_MODEL_ON_NODES If true, material properties are stored on nodes;
 *                           otherwise they are stored per element.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class SEMsolverAcoustoElastic : public Solver
{
 public:
  using AcousticSolverType =
      SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
                enums::physicType::kAcoustic>;
  using ElasticSolverType =
      SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
                enums::physicType::kElastic>;
  using DataType = SEMsolverDataAcoustoElastic;

  SEMsolverAcoustoElastic() = default;
  ~SEMsolverAcoustoElastic() = default;

  // --- Solver interface ---

  int getNumComponents() const override { return 4; }  // p, ux, uy, uz

  /// @brief Returns the acoustic sub-solver mass matrix for DD synchronization.
  VECTOR_REAL_VIEW& getMassMatrixAcoustic() override
  {
    return m_acoustic_solver_.getMassMatrixAcoustic();
  }

  /// @brief Returns the elastic sub-solver mass matrix for DD synchronization.
  VECTOR_REAL_VIEW& getMassMatrixElastic() override
  {
    return m_elastic_solver_.getMassMatrixElastic();
  }

  VECTOR_REAL_VIEW& getDampingMatrix(int c) override
  {
    if (c == 0) return m_acoustic_solver_.getDampingMatrix(0);
    return m_elastic_solver_.getDampingMatrix(c - 1);
  }

  VECTOR_REAL_VIEW& getForceVector(int c) override
  {
    if (c == 0) return m_acoustic_solver_.getForceVector(0);
    return m_elastic_solver_.getForceVector(c - 1);
  }

  void computeFEInit(model::ModelApi<float, int>& mesh,
                     const std::array<float, 3>& sponge_size,
                     const bool surface_sponge,
                     const float taper_delta) override;

  /**
   * @brief Perform one coupled time step (serial / non-distributed mode).
   *
   * Implements the staggered explicit scheme:
   *   elastic forces → elastic update → acoustic forces → acoustic update,
   * with interface coupling applied between the two sub-steps.
   */
  void computeOneStep(const float& dt, const int& timeSample,
                      DataStruct& data) override;

  void computeForces(const float& dt, const int& timeSample,
                     DataStruct& data) override;

  void updateSolution(const float& dt, DataStruct& data) override;

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;
  void computeDampingMatrix() override;

  void outputSolutionValues(const int& t, int& e,
                            const VECTOR_REAL_VIEW& field,
                            const char* fieldName) override;

  void setAnisotropyType(model::AnisotropyType type) override
  {
    m_elastic_solver_.setAnisotropyType(type);
  }

  // --- Accessors for diagnostics ---

  /// Number of acoustic elements detected in the mesh.
  int getNumAcousticElements() const { return num_acoustic_elements_; }

  /// Number of elastic elements detected in the mesh.
  int getNumElasticElements() const { return num_elastic_elements_; }

  /// Number of interface nodes (adjacent to both domains).
  int getNumInterfaceNodes() const { return num_interface_nodes_; }

  // -----------------------------------------------------------------------
  // Methods containing GPU kernels (LOOPHEAD / MAINLOOPHEAD expand to
  // extended __device__ lambdas under CUDA; those must reside in public
  // member functions — CUDA constraint on extended lambdas).
  // -----------------------------------------------------------------------

  /**
   * @brief Identify interface nodes (adjacent to both acoustic and elastic
   * elements).
   */
  void TagNodes();

  /**
   * @brief Compute per-node interface coupling coefficients.
   *
   * For each interface face (between an acoustic and an elastic element),
   * integrates the outward normal over the face using the same quadrature as
   * the damping matrix (INTEGRAL_TYPE::computeDampingTerm).  The result is
   * an area-weighted normal stored in m_coupling_coeff_x/y/z_ at each global
   * node, ready for use in ApplyCouplingAcousticToElastic and
   * ApplyCouplingElasticToAcoustic.
   *
   * The face geometry is queried via faceNormal() and works for any mesh.
   */
  void ComputeInterfaceCouplingCoefficients();

  /**
   * @brief Apply acoustic pressure as traction on elastic interface nodes.
   *
   * Applied POST-Verlet (GEOS style): directly modifies u^{n+1} stored in
   * elastic getPreviousField().  Explicit dt² and M_e factors are visible.
   *
   * @param dt   Time step (used to compute dt²).
   * @param data Coupled solver data containing both wavefields.
   */
  void ApplyCouplingAcousticToElastic(float dt, const DataType& data);

  /**
   * @brief Apply elastic normal acceleration as source on acoustic interface
   * nodes.
   *
   * @param dt   Time step (used to compute acceleration from displacement).
   * @param data Coupled solver data containing both wavefields.
   */
  void ApplyCouplingElasticToAcoustic(float dt, const DataType& data);

 private:
  AcousticSolverType m_acoustic_solver_;  ///< Acoustic sub-solver
  ElasticSolverType m_elastic_solver_;    ///< Elastic sub-solver

  MESH_TYPE m_mesh_;  ///< Local copy of the mesh

  /// Per-element type tag (kElementTypeAcoustic or kElementTypeElastic).
  VECTOR_INT_VIEW m_element_type_;

  /// Boolean mask: true if the node lies on the fluid–solid interface.
  VECTOR_INT_VIEW m_is_interface_node_;

  /// Index map from global node index to interface node index (-1 if not).
  VECTOR_INT_VIEW m_interface_node_index_;

  /// Area-weighted interface coupling coefficient — X component.
  /// Defined for every global node; zero for non-interface nodes.
  VECTOR_REAL_VIEW m_coupling_coeff_x_;
  /// Area-weighted interface coupling coefficient — Y component.
  VECTOR_REAL_VIEW m_coupling_coeff_y_;
  /// Area-weighted interface coupling coefficient — Z component.
  VECTOR_REAL_VIEW m_coupling_coeff_z_;

  /// Elastic displacement at time n-1 (x component), saved before each
  /// elastic Verlet step so that the E→A coupling can use the exact GEOS
  /// finite-difference formula (u^{n+1} - 2*u^n + u^{n-1}) / M_f.
  VECTOR_REAL_VIEW m_ux_nm1_;
  VECTOR_REAL_VIEW m_uy_nm1_;
  VECTOR_REAL_VIEW m_uz_nm1_;

  int num_acoustic_elements_{0};  ///< Count of acoustic elements
  int num_elastic_elements_{0};   ///< Count of elastic elements
  int num_interface_nodes_{0};    ///< Count of interface nodes

  /// Shear-modulus threshold below which an element is classified as acoustic.
  static constexpr float kMuTolerance = 1.0e-6f;

  /**
   * @brief Classify each element as acoustic or elastic based on shear modulus.
   *
   * An element is acoustic (fluid) when mu = rho * vs^2 < kMuTolerance.
   */
  void TagElements();
};

}  // namespace fe
}  // namespace solver

#include "sem_solver_acoustoelastic_impl.h"

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_H_
