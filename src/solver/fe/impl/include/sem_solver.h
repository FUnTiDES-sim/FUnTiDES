#ifndef SRC_SOLVER_FE_API_INCLUDE_SEMSOLVER_H_
#define SRC_SOLVER_FE_API_INCLUDE_SEMSOLVER_H_

#include <array>
#include <cmath>
#include <stdexcept>

#include "data_type.h"
#include "model.h"
#include "parallel_topology.h"
#include "physics_traits.h"
#include "sem_enums.h"
#include "sem_solver_base.h"

namespace solver
{
namespace fe
{

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, enums::physicType PHYSICS>
class SEMsolver : public SEMSolverBase
{
 public:
  using Traits = PhysicsTraits<PHYSICS>;
  using DataType = SEMsolverData<PHYSICS>;

  static constexpr int kNumFields = Traits::kNumFields;
  static constexpr int kNumRhs = Traits::kNumRhsComponents;

  SEMsolver() = default;
  ~SEMsolver() = default;

  // --- Implementation of DD Interface ---

  int getNumComponents() const override { return kNumFields; }

  VECTOR_REAL_VIEW& getMassMatrix() override { return massMatrixGlobal_; }

  VECTOR_REAL_VIEW& getForceVector(int c) override
  {
    return workVectorsGlobal_[c];
  }

  const utils::ParallelTopology& getTopology() const override
  {
    return m_topology;
  }

  // -------------------------------------

  void computeFEInit(model::ModelApi<float, int>& mesh, int rank, int size,
                     float origin_x, float local_lx,
                     const std::array<float, 3>& sponge_size,
                     const bool surface_sponge,
                     const float taper_delta) override;

  // Split-phase methods for DD
  void computeForces(const float& dt, const int& timeSample,
                     DataStruct& data) override;
  void updateSolution(const float& dt, const int& i1, const int& i2,
                      DataStruct& data) override;

  /**
   * @brief Legacy/Serial wrapper.
   * @throws std::runtime_error if called in distributed mode.
   */
  void computeOneStep(const float& dt, const int& timeSample,
                      DataStruct& data) override
  {
    if (m_topology.isDistributed())
    {
      throw std::runtime_error(
          "computeOneStep called in distributed mode. Use computeForces() -> "
          "synchronize() -> updateSolution().");
    }

    auto& myData = dynamic_cast<DataType&>(data);

    computeForces(dt, timeSample, data);
    updateSolution(dt, myData.m_i1, myData.m_i2, data);
  }

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;

  void outputSolutionValues(const int& indexTimeStep, int& i1,
                            int& myElementSource, const ARRAY_REAL_VIEW& field,
                            const char* fieldName) override;

  void applyRHSTerm(int timeSample, float dt, int i2, const DataType& data);
  void computeElementContributions(int i2, const DataType& data);
  void updateFields(float dt, int i1, int i2, DataType& data);

  template <physicType P = PHYSICS,
            typename = std::enable_if_t<P == enums::physicType::kElastic>>
  PROXY_HOST_DEVICE void computeCMatrix(float vp, float vs, float rho,
                                        float delta, float epsilon, float gamma,
                                        float phi, float theta,
                                        float (&C)[6][6]) const;

 private:
  MESH_TYPE m_mesh;
  utils::ParallelTopology m_topology;

  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;

  INTEGRAL_TYPE myQkIntegrals_;

  VECTOR_REAL_VIEW spongeTaperCoeff_;
  VECTOR_REAL_VIEW massMatrixGlobal_;
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

#endif  // SRC_SOLVER_FE_API_INCLUDE_SEMSOLVER_H_
