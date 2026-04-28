#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "data_type.h"
#include "dg_penalty.h"
#include "dg_solver_data.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "parallel_topology.h"
#include "physics_traits.h"
#include "physics_traits_acoustic.h"

#include <typeinfo>


#include "sem_enums.h"
#include "solver.h"

namespace solver {
namespace fe {

using physicType = utils::enums::physicType;

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
class DGsolver : public Solver {
 public:

  using Traits = PhysicsTraits<PHYSICS>;

  using DataType = DGsolverDataAcoustic;

  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;

  DGsolver() = default;
  ~DGsolver() = default;

  // --- Implementation of DD Interface ---

  int getNumComponents() const override { return kNumFields; }

  // --- Mandatory overrides for Solver interface ---

  void initFEarrays() override {
    // TODO: Implement initialization logic
  }

  void allocateFEarrays() override {
    // TODO: Implement allocation logic
  }

  void initSpongeValues() override {
    // TODO: Implement sponge layer initialization
  }

  void resetGlobalVectors(int numNodes) override {
    // TODO: Implement reset logic
  }

  void computeGlobalMassMatrix() override {
    // TODO: Implement mass matrix calculation
  }

  void computeDampingMatrix() override {
    // TODO: Implement damping matrix calculation
  }

  void outputSolutionValues(const int& t, int& e, const VECTOR_REAL_VIEW& field, const char* fieldName) override {}
     
  void outputSolutionValues(const int& t, int& e, const ARRAY_REAL_VIEW& field, const char* fieldName) override;

  VECTOR_REAL_VIEW& getMassMatrixAcoustic() override {
    throw std::runtime_error("getMassMatrixAcoustic not implemented for DG");
  }

  VECTOR_REAL_VIEW& getMassMatrixElastic() override {
    throw std::runtime_error("getMassMatrixElastic not implemented for DG");
  }

  VECTOR_REAL_VIEW& getDampingMatrix(int c) override {
    throw std::runtime_error("getDampingMatrix not implemented for DG");
  }

  VECTOR_REAL_VIEW& getForceVector(int component) override {
    throw std::runtime_error("getForceVector not implemented for DG");
  }

  void setAnisotropyType(model::AnisotropyType type) override {
    // TODO: Implement anisotropy setting
  }

  void setSLSAttenuation(const VECTOR_REAL_VIEW& reference_frequencies,
                         const VECTOR_REAL_VIEW& anelasticity_coefficients = VECTOR_REAL_VIEW()) override {
    // TODO: Implement SLS attenuation setting
  }

  // --- Core solver methods ---

  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override;

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

  /**
   * @brief Apply external forcing to the global fields.
   */
  void applyRHSTerm(int timeSample, float dt, const DataType& data);

  /**
   * @brief Update the global solution fields at interior nodes.
   */
  void updateFields(float dt, const DataType& data);

 private:
  MESH_TYPE m_mesh;
  model::FaceConnectivityUnstruct<float, int> m_face_connectivity_;
  real_t m_penalty_factor_ = 10.0f;

  std::array<VECTOR_REAL_VIEW, kNumFields> rhsTermGlobal;

  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  static constexpr int knumNodesPerFace = (ORDER + 1) * (ORDER + 1);
};

// Backward Compatibility Aliases
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
using DGsolverAcoustic =
    DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, utils::enums::physicType::kAcoustic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_H_