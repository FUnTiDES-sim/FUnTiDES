#ifndef SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#define SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
#include <array>
#include <cmath>
#include <stdexcept>

#include "data_type.h"
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

  VECTOR_REAL_VIEW& getMassMatrix() override { return massMatrixGlobal_; }

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
   * @brief Update the global solution fields at interior nodes.
   *
   * @param dt Delta time for this iteration.
   * @param data Data structure containing solution fields.
   */
  void updateFields(float dt, const DataType& data);

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

 private:
  MESH_TYPE m_mesh;

  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;

  INTEGRAL_TYPE myQkIntegrals_;

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
#endif  // SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_H_
