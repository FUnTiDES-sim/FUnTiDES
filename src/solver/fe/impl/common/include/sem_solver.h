#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "data_type.h"
#include "face_connectivity_unstruct.h"
#include "mesh_type_traits.h"
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
class SEMsolver : public Solver {
 public:
  using Traits = PhysicsTraits<PHYSICS>;
  using DataType = SEMsolverData<PHYSICS>;

  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;

  SEMsolver() = default;
  ~SEMsolver() = default;

  // --- Implementation of DD Interface ---

  int getNumComponents() const override { return kNumFields; }

  VECTOR_REAL_VIEW& getMassMatrixAcoustic() override { return massMatrixGlobal_; }
  VECTOR_REAL_VIEW& getMassMatrixElastic() override { return massMatrixGlobal_; }

  VECTOR_REAL_VIEW& getDampingMatrix(int c) override { return dampingMatrixGlobal_[c]; }

  VECTOR_REAL_VIEW& getForceVector(int c) override { return workVectorsGlobal_[c]; }

  // -------------------------------------

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

  // Split-phase methods for DD
  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override;
  void updateSolution(const float& dt, DataStruct& data) override;

  /**
   * @brief Legacy/Serial wrapper.
   * @throws std::runtime_error if called in distributed mode.
   */
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override {
    auto& myData = dynamic_cast<DataType&>(data);
    if (myData.isDistributed) {
      throw std::runtime_error(
          "computeOneStep called in distributed mode. Use computeForces() -> "
          "synchronize() -> updateSolution().");
    }
    computeForces(dt, timeSample, data);
    updateSolution(dt, data);
  }

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;
  void computeDampingMatrix() override;

  /**
   * @brief Assemble mass matrix restricted to elements matching a mask.
   *
   * Identical to computeGlobalMassMatrix() but skips elements where
   * @p elem_mask[e] != @p active_value.  Used by SEMsolverAcoustoElastic to
   * delegate domain-specific mass matrix assembly to each sub-solver without
   * duplicating the quadrature kernel.
   *
   * @param elem_mask   Per-element integer tag array (size = nElements).
   * @param active_value Only elements with this tag value are accumulated.
   */
  void computeGlobalMassMatrixMasked(const VECTOR_INT_VIEW& elem_mask, int active_value);

  /**
   * @brief Assemble damping matrix restricted to elements matching a mask.
   *
   * Parallel variant of computeDampingMatrix() that skips elements where
   * @p elem_mask[e] != @p active_value.  Used by SEMsolverAcoustoElastic to
   * delegate domain-specific assembly to each sub-solver without duplicating
   * the quadrature kernel.
   *
   * @param elem_mask   Per-element integer tag array (size = nElements).
   * @param active_value Only elements with this tag value are accumulated.
   */
  void computeDampingMatrixMasked(const VECTOR_INT_VIEW& elem_mask, int active_value);

  void outputSolutionValues(const int& t, int& e, const VECTOR_REAL_VIEW& field, const char* fieldName) override;

  /**
   * @brief Apply external forcing to the global fields.
   *
   * @param timeSample Current time sample index.
   * @param dt Delta time for this iteration.
   * @param data Data structure containing RHS terms and fields.
   */
  void applyRHSTerm(int timeSample, float dt, const DataType& data);

  /**
   * @brief Assemble local element contributions to global FE vectors.
   *
   * @param data Data structure containing solution fields.
   */
  void computeElementContributions(const DataType& data);

  /**
   * @brief Assemble local element contributions restricted to elements
   *        matching a mask.
   *
   * Identical to computeElementContributions() but skips elements where
   * @p elem_mask[e] != @p active_value.  Used by SEMsolverAcoustoElastic to
   * delegate domain-specific stiffness assembly to each sub-solver without
   * duplicating the quadrature kernel.
   *
   * @param data         Data structure containing solution fields.
   * @param elem_mask    Per-element integer tag array (size = nElements).
   * @param active_value Only elements with this tag value are accumulated.
   */
  void computeElementContributionsMasked(const DataType& data, const VECTOR_INT_VIEW& elem_mask, int active_value);

  /**
   * @brief Assemble stiffness into workVectorsGlobal_ for a compact list of
   *        elements, skipping all others.
   * @param data      Data structure containing solution fields.
   * @param elem_list Compact array of element indices to process.
   * @param n_elems   Number of entries in @p elem_list.
   */
  void computeElementContributionsFromList(const DataType& data, const VECTOR_INT_VIEW& elem_list, int n_elems);

  /**
   * @brief Update the global solution fields at interior nodes.
   *
   * @param dt Delta time for this iteration.
   * @param data Data structure containing solution fields.
   */
  void updateFields(float dt, const DataType& data);

  /**
   * @brief Run Verlet update only for a compact subset of nodes.
   * @param dt        Time step.
   * @param data      Wavefield data.
   * @param node_list Compact array of node indices to update.
   * @param n_nodes   Number of entries in @p node_list.
   */
  void updateFieldsFromList(float dt, const DataType& data, const VECTOR_INT_VIEW& node_list, int n_nodes);

  /**
   * @brief Read-only access to the f-th work (force) vector.
   *
   * After computeElementContributions[FromList], workVectorsGlobal_[f] holds
   * the assembled RHS + K*u^n for component f.  This accessor lets the
   * AcoustoElastic solver read the elastic force without friendship.
   *
   * @param f Component index (0 <= f < kNumFields).
   * @return  Const reference to the view.
   */
  const VECTOR_REAL_VIEW& getForceVector(int f) const { return workVectorsGlobal_[f]; }

  void computeElementContributions_Acoustic(const DataType& data);
  void computeElementContributions_Iso(const DataType& data);
  void computeElementContributions_Vti(const DataType& data);
  void computeElementContributions_Tti(const DataType& data);

  void computeAttenuationContributions(const DataType& data);
  void computeAttenuationContributionsAcoustic(const DataType& data);
  void computeAttenuationContributionsElastic(const DataType& data);

  /**
   * @brief Compute the elasticity matrix at a given node (elastic only).
   *
   * This method is only compiled for elastic physics.
   *
   * @param vp P-wave velocity.
   * @param vs S-wave velocity.
   * @param rho Density.
   * @param delta Thomsen parameter delta.
   * @param epsilon Thomsen parameter epsilon.
   * @param gamma Thomsen parameter gamma.
   * @param phi Azimuthal angle (radians).
   * @param theta Dip angle (radians).
   * @param C Output 6x6 elasticity matrix.
   */
  template <physicType P = PHYSICS, typename = std::enable_if_t<P == utils::enums::physicType::kElastic>>
  static PROXY_HOST_DEVICE void computeCMatrix(float vp, float vs, float rho, float delta, float epsilon, float gamma,
                                               float phi, float theta, float (&C)[6][6]);

  /**
   * @brief Set the anisotropy type for the solver.
   */
  void setAnisotropyType(model::AnisotropyType type) { anisotropyType_ = type; }

  void setSLSAttenuation(const VECTOR_REAL_VIEW& reference_frequencies,
                         const VECTOR_REAL_VIEW& anelasticity_coefficients = VECTOR_REAL_VIEW()) override {
    attenuationEnabled_ = reference_frequencies.extent(0) > 0;
    nSls_ = static_cast<int>(reference_frequencies.extent(0));
    if (!attenuationEnabled_) {
      nSls_ = 0;
      slsReferenceAngularFrequencies_ = VECTOR_REAL_VIEW();
      slsAnelasticityCoefficients_ = VECTOR_REAL_VIEW();
      return;
    }

    slsReferenceAngularFrequencies_ = allocateVector<VECTOR_REAL_VIEW>(nSls_, "slsReferenceAngularFrequencies");
    for (int i = 0; i < nSls_; ++i) {
      slsReferenceAngularFrequencies_[i] = reference_frequencies[i];
    }

    slsAnelasticityCoefficients_ = allocateVector<VECTOR_REAL_VIEW>(nSls_, "slsAnelasticityCoefficients");
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

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;
  model::AnisotropyType anisotropyType_;

  INTEGRAL_TYPE myQkIntegrals_;

  // Mask state used by computeElementContributionsMasked.
  bool m_mask_enabled_ = false;
  VECTOR_INT_VIEW m_element_mask_;
  int m_mask_active_value_ = 0;

  // List state used by computeElementContributionsFromList.
  bool m_list_mode_ = false;
  VECTOR_INT_VIEW m_elem_list_;
  int m_n_elem_list_ = 0;

  // Node list state used by updateFieldsFromList.
  bool m_node_list_mode_ = false;
  VECTOR_INT_VIEW m_node_list_;
  int m_n_node_list_ = 0;

  VECTOR_REAL_VIEW spongeTaperCoeff_;
  VECTOR_REAL_VIEW massMatrixGlobal_;
  std::array<VECTOR_REAL_VIEW, kNumFields> dampingMatrixGlobal_;
  std::array<VECTOR_REAL_VIEW, kNumFields> workVectorsGlobal_;

  typename INTEGRAL_TYPE::Metrics3D m_struct_metrics_{};
  bool m_has_struct_metrics_ = false;

  bool attenuationEnabled_ = false;
  int nSls_ = 0;
  VECTOR_REAL_VIEW slsReferenceAngularFrequencies_;
  VECTOR_REAL_VIEW slsAnelasticityCoefficients_;
  std::array<VECTOR_REAL_VIEW, kNumFields> attenuationWorkVectorsGlobal_;
  std::array<ARRAY_REAL_VIEW, kNumFields> attenuationMemoryVariables_;
};

// Backward Compatibility Aliases
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using SEMsolverAcoustic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kAcoustic>;

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using SEMsolverElastic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kElastic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
