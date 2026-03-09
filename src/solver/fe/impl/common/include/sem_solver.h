#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#include <array>
#include <cmath>
#include <stdexcept>

#include "data_type.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "parallel_topology.h"
#include "physics_traits.h"
#include "sem_enums.h"
#include "sem_solver_data.h"
#include "solver.h"

namespace solver
{
namespace fe
{

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, enums::physicType PHYSICS>
class SEMsolver : public Solver
{
 public:
  using Traits = PhysicsTraits<PHYSICS>;
  using DataType = SEMsolverData<PHYSICS>;

  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;

  SEMsolver() = default;
  ~SEMsolver() = default;

  // --- Implementation of DD Interface ---

  int getNumComponents() const override { return kNumFields; }

  VECTOR_REAL_VIEW& getMassMatrixAcoustic() override
  {
    return massMatrixGlobal_;
  }
  VECTOR_REAL_VIEW& getMassMatrixElastic() override
  {
    return massMatrixGlobal_;
  }

  VECTOR_REAL_VIEW& getDampingMatrix(int c) override
  {
    return dampingMatrixGlobal_[c];
  }

  VECTOR_REAL_VIEW& getForceVector(int c) override
  {
    return workVectorsGlobal_[c];
  }

  // -------------------------------------

  void computeFEInit(model::ModelApi<float, int>& mesh,
                     const std::array<float, 3>& sponge_size,
                     const bool surface_sponge,
                     const float taper_delta) override;

  // Split-phase methods for DD
  void computeForces(const float& dt, const int& timeSample,
                     DataStruct& data) override;
  void updateSolution(const float& dt, DataStruct& data) override;

  /**
   * @brief Legacy/Serial wrapper.
   * @throws std::runtime_error if called in distributed mode.
   */
  void computeOneStep(const float& dt, const int& timeSample,
                      DataStruct& data) override
  {
    auto& myData = dynamic_cast<DataType&>(data);
    if (myData.isDistributed)
    {
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
  void computeGlobalMassMatrixMasked(const VECTOR_INT_VIEW& elem_mask,
                                     int active_value);

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
  void computeDampingMatrixMasked(const VECTOR_INT_VIEW& elem_mask,
                                  int active_value);

  void outputSolutionValues(const int& t, int& e, const VECTOR_REAL_VIEW& field,
                            const char* fieldName) override;

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
  void computeElementContributionsMasked(const DataType& data,
                                         const VECTOR_INT_VIEW& elem_mask,
                                         int active_value);

  /**
   * @brief Assemble stiffness into workVectorsGlobal_ for a compact list of
   *        elements, skipping all others.
   * @param data      Data structure containing solution fields.
   * @param elem_list Compact array of element indices to process.
   * @param n_elems   Number of entries in @p elem_list.
   */
  void computeElementContributionsFromList(const DataType& data,
                                           const VECTOR_INT_VIEW& elem_list,
                                           int n_elems);

  /**
   * @brief Update the global solution fields at interior nodes.
   *
   * @param dt Delta time for this iteration.
   * @param data Data structure containing solution fields.
   */
  void updateFields(float dt, const DataType& data);

  void computeElementContributions_Iso(const DataType& data);
  void computeElementContributions_Vti(const DataType& data);
  void computeElementContributions_Tti(const DataType& data);

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
  template <physicType P = PHYSICS,
            typename = std::enable_if_t<P == enums::physicType::kElastic>>
  PROXY_HOST_DEVICE void computeCMatrix(float vp, float vs, float rho,
                                        float delta, float epsilon, float gamma,
                                        float phi, float theta,
                                        float (&C)[6][6]) const;

  /**
   * @brief Set the anisotropy type for the solver.
   */
  void setAnisotropyType(model::AnisotropyType type) { anisotropyType_ = type; }

 private:
  MESH_TYPE m_mesh;

  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

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

  VECTOR_REAL_VIEW spongeTaperCoeff_;
  VECTOR_REAL_VIEW massMatrixGlobal_;
  std::array<VECTOR_REAL_VIEW, kNumFields> dampingMatrixGlobal_;
  std::array<VECTOR_REAL_VIEW, kNumFields> workVectorsGlobal_;
};

// Backward Compatibility Aliases
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
using SEMsolverAcoustic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
              enums::physicType::kAcoustic>;

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
using SEMsolverElastic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
              enums::physicType::kElastic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
