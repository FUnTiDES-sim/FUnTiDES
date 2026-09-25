#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "data_type.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "parallel_topology.h"
#include "physics_traits.h"
#include "physics_traits_acoustic.h"
#include "physics_traits_elastic.h"
#include "sem_enums.h"
#include "sem_solver_data.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Spectral-element solver advancing one physics with an explicit Verlet time scheme.
 *
 * Owns the global mass, damping and work vectors indexed by global node. Each time step
 * assembles the element contributions into the work vectors, then updates the wavefield.
 * The step can be split (computeForces, then updateSolutionForward) so that a caller can
 * synchronize shared nodes in between.
 *
 * @tparam ORDER Polynomial order of the elements.
 * @tparam INTEGRAL_TYPE Element integral back-end.
 * @tparam MESH_TYPE Mesh type held by the solver.
 * @tparam IS_MODEL_ON_NODES True if the model properties are stored per node, false if per element.
 * @tparam PHYSICS Physics solved.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
class SEMsolver : public Solver {
 public:
  using Traits = PhysicsTraits<PHYSICS>;
  using DataType = SEMsolverData<PHYSICS>;

  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;  ///< Wavefield components.
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;    ///< Source components.

  SEMsolver() = default;
  ~SEMsolver() = default;

  int getNumComponents() const override { return kNumFields; }

  vectorReal& getMassMatrixAcoustic() override { return massMatrixGlobal_; }
  vectorReal& getMassMatrixElastic() override { return massMatrixGlobal_; }

  vectorReal& getDampingMatrix(int c) override { return dampingMatrixGlobal_[c]; }

  /// @return Sponge taper coefficients.
  /// @todo VERIFY: size and indexing of the taper vector (per node?).
  vectorReal& getSpongeTaperCoeff() { return spongeTaperCoeff_; }

  vectorReal& getForceVector(int c) override { return workVectorsGlobal_[c]; }

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override;
  void updateSolutionForward(const float& dt, DataStruct& data) override;
  void updateSolutionBackward(const float& dt, DataStruct& data) override;

  /**
   * @brief Advance one forward time step without any inter-rank synchronization.
   *
   * Calls computeForces() then updateSolutionForward().
   * @throws std::runtime_error if the data is flagged as distributed.
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

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;
  void computeDampingMatrix() override;

  /**
   * @brief Assemble the mass matrix over the elements carrying a given tag.
   *
   * Same as computeGlobalMassMatrix() but skips elements whose tag differs from @p active_value.
   *
   * @param elem_mask Per-element integer tag, size numElements.
   * @param active_value Only elements with this tag are accumulated.
   */
  void computeGlobalMassMatrixMasked(const vectorInt& elem_mask, int active_value);

  /**
   * @brief Assemble the damping matrix over the elements carrying a given tag.
   *
   * Same as computeDampingMatrix() but skips elements whose tag differs from @p active_value.
   *
   * @param elem_mask Per-element integer tag, size numElements.
   * @param active_value Only elements with this tag are accumulated.
   */
  void computeDampingMatrixMasked(const vectorInt& elem_mask, int active_value);

  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override;

  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override {};

  /**
   * @brief Add the external source term to the global work vectors.
   *
   * @param timeSample Current time sample index.
   * @param dt Time step.
   * @param data Solver data holding the source and the fields.
   */
  void applyRHSTerm(int timeSample, float dt, const DataType& data);

  /**
   * @brief Assemble the element stiffness contributions into the global work vectors.
   *
   * @param data Solver data holding the wavefield.
   */
  void computeElementContributions(const DataType& data);

  /**
   * @brief Assemble the element stiffness contributions over the elements carrying a given tag.
   *
   * Same as computeElementContributions() but skips elements whose tag differs from @p active_value.
   *
   * @param data Solver data holding the wavefield.
   * @param elem_mask Per-element integer tag, size numElements.
   * @param active_value Only elements with this tag are accumulated.
   */
  void computeElementContributionsMasked(const DataType& data, const vectorInt& elem_mask, int active_value);

  /**
   * @brief Assemble the element stiffness contributions over an explicit list of elements.
   *
   * @param data Solver data holding the wavefield.
   * @param elem_list Compact array of element indices to process.
   * @param n_elems Number of valid entries in @p elem_list.
   */
  void computeElementContributionsFromList(const DataType& data, const vectorInt& elem_list, int n_elems);

  /**
   * @brief Verlet update of the wavefield at all nodes, forward in time.
   *
   * @param dt Time step.
   * @param data Solver data holding the wavefield.
   */
  void updateFieldsForward(float dt, const DataType& data);

  /**
   * @brief Verlet update of the wavefield at all nodes for adjoint time stepping.
   *
   * Same as updateFieldsForward() but writes the result to the prevprev buffer instead of prev.
   *
   * @param dt Time step.
   * @param data Solver data holding the wavefield; the prevprev buffer must be allocated.
   */
  void updateFieldsBackward(float dt, const DataType& data);

  /**
   * @brief Forward Verlet update restricted to a list of nodes.
   *
   * @param dt Time step.
   * @param data Solver data holding the wavefield.
   * @param node_list Compact array of node indices to update.
   * @param n_nodes Number of valid entries in @p node_list.
   */
  void updateFieldsFromListForward(float dt, const DataType& data, const vectorInt& node_list, int n_nodes);

  /**
   * @brief Backward Verlet update restricted to a list of nodes.
   *
   * @param dt Time step.
   * @param data Solver data holding the wavefield; the prevprev buffer must be allocated.
   * @param node_list Compact array of node indices to update.
   * @param n_nodes Number of valid entries in @p node_list.
   */
  void updateFieldsFromListBackward(float dt, const DataType& data, const vectorInt& node_list, int n_nodes);

  /**
   * @brief Read-only access to the f-th work vector.
   *
   * After the element contributions are assembled, it holds the source plus the stiffness
   * term applied to the current field, for component f.
   *
   * @param f Component index, 0 <= f < kNumFields.
   * @return Const reference to the view.
   */
  const vectorReal& getForceVector(int f) const { return workVectorsGlobal_[f]; }

  /// Acoustic stiffness assembly: default, flat (one thread per element) and GEMM variants.
  void computeElementContributions_Acoustic(const DataType& data);
  void computeElementContributions_Acoustic_Flat(const DataType& data);
  void computeElementContributions_Acoustic_Gemm(const DataType& data);

  /**
   * @brief Highest order still served by the one-thread-per-element kernels.
   *
   * The team kernels give a whole warp to the kPointsPerElement quadrature points of an
   * element, so at order 1 (8 points) most of the warp idles at the barriers while the flat
   * kernel keeps every thread busy on its own element. The crossover is hardware-dependent;
   * re-measure it when moving to another GPU.
   */
  static constexpr int kMaxOrderForFlatElastic = 1;

  /// Elastic stiffness assembly per anisotropy type (Iso, Vti, Tti): default, flat and team variants.
  void computeElementContributions_Iso(const DataType& data);
  void computeElementContributions_Iso_Flat(const DataType& data);
  void computeElementContributions_Iso_Team(const DataType& data);
  void computeElementContributions_Vti(const DataType& data);
  void computeElementContributions_Vti_Flat(const DataType& data);
  void computeElementContributions_Vti_Team(const DataType& data);
  void computeElementContributions_Tti(const DataType& data);
  void computeElementContributions_Tti_Flat(const DataType& data);
  void computeElementContributions_Tti_Team(const DataType& data);

  /// Add the attenuation (SLS) contributions to the work vectors; no effect unless attenuation is enabled.
  void computeAttenuationContributions(const DataType& data);
  void computeAttenuationContributionsAcoustic(const DataType& data);
  void computeAttenuationContributionsElastic(const DataType& data);

  /**
   * @brief Compute the 6x6 elasticity matrix at a node. Elastic physics only.
   *
   * @param vp P-wave velocity.
   * @param vs S-wave velocity.
   * @param rho Density.
   * @param delta Thomsen parameter delta.
   * @param epsilon Thomsen parameter epsilon.
   * @param gamma Thomsen parameter gamma.
   * @param phi Azimuthal angle (radians).
   * @param theta Dip angle (radians).
   * @param[out] C 6x6 elasticity matrix.
   * @todo VERIFY: Voigt component order of C.
   */
  template <physicType P = PHYSICS, typename = std::enable_if_t<P == utils::enums::physicType::kElastic>>
  static PROXY_HOST_DEVICE void computeCMatrix(float vp, float vs, float rho, float delta, float epsilon, float gamma,
                                               float phi, float theta, float (&C)[6][6]);

  /**
   * @brief Build the per-node TTI elasticity tensors, once per model.
   *
   * The model is constant during the time loop, so the rotated tensor of a node is the same at
   * every step. Building it once keeps the rotation out of the stiffness kernel. No-op unless
   * the physics is elastic and the model lives on nodes; cheap to call repeatedly.
   */
  void precomputeTtiTensorsOnNodes();

  /// Select the anisotropy type used by the elastic stiffness kernels.
  void setAnisotropyType(model::AnisotropyType type) { anisotropyType_ = type; }

  /**
   * @brief Enable or disable SLS attenuation.
   *
   * An empty @p reference_frequencies disables attenuation.
   *
   * @param reference_frequencies Reference angular frequencies of the SLS mechanisms.
   * @param anelasticity_coefficients One coefficient per mechanism; if empty, every coefficient
   *        is set to -1.
   * @throws std::runtime_error if the two arrays have different sizes.
   * @todo VERIFY: meaning of the default coefficient -1, and whether the frequencies are angular (rad/s) or in Hz.
   */
  void setSLSAttenuation(const vectorReal& reference_frequencies,
                         const vectorReal& anelasticity_coefficients = vectorReal()) override {
    attenuationEnabled_ = reference_frequencies.extent(0) > 0;
    nSls_ = static_cast<int>(reference_frequencies.extent(0));
    if (!attenuationEnabled_) {
      nSls_ = 0;
      slsReferenceAngularFrequencies_ = vectorReal();
      slsAnelasticityCoefficients_ = vectorReal();
      return;
    }

    slsReferenceAngularFrequencies_ = allocateVector<vectorReal>(nSls_, "slsReferenceAngularFrequencies");
    for (int i = 0; i < nSls_; ++i) {
      slsReferenceAngularFrequencies_[i] = reference_frequencies[i];
    }

    slsAnelasticityCoefficients_ = allocateVector<vectorReal>(nSls_, "slsAnelasticityCoefficients");
    if (anelasticity_coefficients.extent(0) == 0) {
      for (int i = 0; i < nSls_; ++i) {
        slsAnelasticityCoefficients_[i] = -1.0f;
      }
    } else {
      if (static_cast<int>(anelasticity_coefficients.extent(0)) != nSls_) {
        throw std::runtime_error(
            "SLS anelasticity coefficients must match reference frequencies "
            "size");
      }
      for (int i = 0; i < nSls_; ++i) {
        slsAnelasticityCoefficients_[i] = anelasticity_coefficients[i];
      }
    }
  }

 private:
  MESH_TYPE m_mesh;

  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  /// The rotated tensor is symmetric, so only its upper triangle (21 entries) is stored.
  static constexpr int kCttiPackedSize = 21;

  vectorReal gemmMetrics_;
  bool gemmMetricsReady_ = false;

  arrayReal cttiNodes_;  ///< Packed rotated TTI tensor per node, kCttiPackedSize entries each.
  bool cttiNodesReady_ = false;

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;
  model::AnisotropyType anisotropyType_;

  INTEGRAL_TYPE myQkIntegrals_;

  // Mask state used by computeElementContributionsMasked.
  bool m_mask_enabled_ = false;
  vectorInt m_element_mask_;
  int m_mask_active_value_ = 0;

  // List state used by computeElementContributionsFromList.
  bool m_list_mode_ = false;
  vectorInt m_elem_list_;
  int m_n_elem_list_ = 0;

  // Node list state used by updateFieldsFromList*.
  bool m_node_list_mode_ = false;
  vectorInt m_node_list_;
  int m_n_node_list_ = 0;

  vectorReal spongeTaperCoeff_;
  vectorReal massMatrixGlobal_;                             ///< Size numNodes.
  std::array<vectorReal, kNumFields> dampingMatrixGlobal_;  ///< One vector of size numNodes per component.
  std::array<vectorReal, kNumFields> workVectorsGlobal_;    ///< One vector of size numNodes per component.

  bool attenuationEnabled_ = false;
  int nSls_ = 0;
  vectorReal slsReferenceAngularFrequencies_;
  vectorReal slsAnelasticityCoefficients_;
  std::array<vectorReal, kNumFields> attenuationWorkVectorsGlobal_;
  std::array<arrayReal, kNumFields> attenuationMemoryVariables_;
};

/// Acoustic specialization of SEMsolver.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using SEMsolverAcoustic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kAcoustic>;

/// Elastic specialization of SEMsolver.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using SEMsolverElastic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kElastic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
