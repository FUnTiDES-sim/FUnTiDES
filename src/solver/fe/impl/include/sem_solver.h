#ifndef SRC_SOLVER_FE_API_INCLUDE_SEMSOLVER_H_
#define SRC_SOLVER_FE_API_INCLUDE_SEMSOLVER_H_

#include <array>
#include <cmath>

#include "data_type.h"
#include "model.h"
#include "physics_traits.h"
#include "sem_enums.h"
#include "sem_solver_base.h"

namespace solver
{
namespace fe
{

//============================================================================
// Unified SEM Solver Class
//============================================================================

/**
 * @brief Unified Spectral Element Method solver for wave propagation.
 *
 * This class implements both acoustic and elastic wave propagation using
 * compile-time physics selection. The physics-specific behavior is handled
 * via `if constexpr` branches, resulting in zero runtime overhead.
 *
 * @tparam ORDER Polynomial order of the spectral elements
 * @tparam INTEGRAL_TYPE Type for numerical integration (basis functions,
 *                       quadrature)
 * @tparam MESH_TYPE Type of the computational mesh
 * @tparam IS_MODEL_ON_NODES True if model parameters are on nodes, false if on
 *                           elements
 * @tparam PHYSICS Physics type (kAcoustic or kElastic)
 */
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

  /**
   * @brief Initialize all finite element structures.
   *
   * Sets up basis functions, integrals, global arrays, and sponge boundaries.
   *
   * @param mesh BaseMesh structure containing the domain information.
   * @param sponge_size Thickness (in elements) of absorbing sponge layers
   *                    in each direction [x, y, z].
   * @param surface_sponge Enable sponge at free surface.
   * @param taper_delta Attenuation parameter for sponge layers.
   */
  void computeFEInit(model::ModelApi<float, int>& mesh,
                     const std::array<float, 3>& sponge_size,
                     const bool surface_sponge,
                     const float taper_delta) override;

  /**
   * @brief Compute one time step of the wave equation solver.
   *
   * @param dt Delta time for this iteration.
   * @param timeSample Current time index into the RHS (source) term.
   * @param data DataStruct containing all necessary arrays.
   */
  void computeOneStep(const float& dt, const int& timeSample,
                      DataStruct& data) override;

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;

  void outputSolutionValues(const int& indexTimeStep, int& i1,
                            int& myElementSource, const ARRAY_REAL_VIEW& field,
                            const char* fieldName) override;

  /**
   * @brief Apply external forcing to the global fields.
   *
   * @param timeSample Current time sample index.
   * @param dt Delta time for this iteration.
   * @param i2 Current field index.
   * @param data Data structure containing RHS terms and fields.
   */
  void applyRHSTerm(int timeSample, float dt, int i2, const DataType& data);

  /**
   * @brief Assemble local element contributions to global FE vectors.
   *
   * @param i2 Current field index.
   * @param data Data structure containing solution fields.
   */
  void computeElementContributions(int i2, const DataType& data);

  /**
   * @brief Update the global solution fields at interior nodes.
   *
   * @param dt Delta time for this iteration.
   * @param i1 Previous time step index.
   * @param i2 Current time step index.
   * @param data Data structure containing solution fields.
   */
  void updateFields(float dt, int i1, int i2, const DataType& data);

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
  MESH_TYPE m_mesh;  ///< Computational mesh

  /// Number of nodes per element
  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float sponge_size_[3];  ///< Sponge layer thickness [x, y, z]
  bool surface_sponge_;   ///< Enable sponge at free surface
  float taper_delta_;     ///< Attenuation parameter

  /// Basis functions and integral objects
  INTEGRAL_TYPE myQkIntegrals_;

  VECTOR_REAL_VIEW spongeTaperCoeff_;  ///< Sponge tapering coefficients
  VECTOR_REAL_VIEW massMatrixGlobal_;  ///< Global mass matrix

  /// Work vectors for accumulating element contributions
  /// Size is kNumFields (1 for acoustic, 3 for elastic)
  std::array<VECTOR_REAL_VIEW, kNumFields> workVectorsGlobal_;
};

//============================================================================
// Backward Compatibility Type Aliases for Solvers
//============================================================================

/**
 * @brief Type alias for acoustic solver (backward compatibility).
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
using SEMsolverAcoustic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
              enums::physicType::kAcoustic>;

/**
 * @brief Type alias for elastic solver (backward compatibility).
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
using SEMsolverElastic =
    SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
              enums::physicType::kElastic>;

}  // namespace fe
}  // namespace solver

#endif  // SRC_SOLVER_FE_API_INCLUDE_SEMSOLVER_H_
