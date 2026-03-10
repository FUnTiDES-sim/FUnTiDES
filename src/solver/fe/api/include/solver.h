#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_SOLVER_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_SOLVER_H_
#include <array>
#include <cmath>
#include <vector>

#include "model.h"
#include "parallel_topology.h"

namespace solver
{
namespace fe
{
/**
 * @brief Base FE solver.
 */
class Solver
{
 public:
  Solver() = default;
  virtual ~Solver() = default;

  struct DataStruct
  {
    // Base structure for data passed to the solver
    virtual ~DataStruct() = default;

    virtual void print() const = 0;
  };

  // Pure virtual function to compute one step of the solver
  virtual void computeOneStep(const float& dt, const int& timeSample,
                              DataStruct& data) = 0;

  /**
   * @brief Initialize all finite element structures:
   * basis functions, integrals, global arrays, and sponge boundaries.
   *
   * @param mesh BaseMesh structure containing the domain information.
   * @param origin_x Local subdomain origin in X (for topology discovery).
   * @param local_lx Local subdomain width in X (for topology discovery).
   * @param sponge_size Thickness (in elements) of absorbing sponge layers
   *                    in each direction [x, y, z] to prevent reflections.
   * @param surface_sponge Enable sponge at free surface (typically false
   *                       for geophysics to preserve natural reflections).
   * @param taper_delta_ Attenuation parameter for sponge layers.
   */
  virtual void computeFEInit(model::ModelApi<float, int>& mesh,
                             const std::array<float, 3>& sponge_size,
                             const bool surface_sponge,
                             const float taper_delta_) = 0;

  /**
   * @brief Initialize arrays required by the finite element solver.
   */
  virtual void initFEarrays() = 0;

  /**
   * @brief Allocate memory for FE-related arrays (mass, stiffness, etc.).
   */
  virtual void allocateFEarrays() = 0;

  /**
   * @brief Initialize sponge (absorbing layer) coefficients.
   */
  virtual void initSpongeValues() = 0;

  /**
   * @brief Reset global FE vectors (mass, stiffness) before accumulation.
   *
   * @param numNodes Total number of global nodes.
   */
  virtual void resetGlobalVectors(int numNodes) = 0;

  /**
   * @brief Compute the global mass matrix.
   *  once at the beginning of the simulation.
   */
  virtual void computeGlobalMassMatrix() = 0;

  /**
   * @brief Compute the global damping matrix.
   *  once at the beginning of the simulation.
   */
  virtual void computeDampingMatrix() = 0;

  /**
   * @brief Outputs solution field values at a specific time step
   *
   * This pure virtual function is responsible for writing or displaying
   * solution values for a given field at a particular time step. It must be
   * implemented by derived classes to define the specific output behavior.
   *
   * @param[in] t The index of the current time step for which
   *              solution values are being output
   * @param[in] e Index of the source element being processed
   * @param[in] field Constant view of the array containing the field values
   *                  to be output
   * @param[in] fieldName Name/identifier of the field being output (as a
   *                      C-string)
   */
  virtual void outputSolutionValues(const int& t, int& e,
                                    const VECTOR_REAL_VIEW& field,
                                    const char* fieldName) = 0;

  // --- Domain Decomposition Interface ---

  /**
   * @brief Get the number of physical field components (e.g., 1 for Acoustic, 3
   * for Elastic).
   */
  virtual int getNumComponents() const = 0;

  /**
   * @brief Access the Global Mass Matrix.
   * Used by the orchestrator to synchronize mass values at boundaries after
   * initialization.
   */
  virtual VECTOR_REAL_VIEW& getMassMatrix() = 0;

  /**
   * @brief Access the Global Damping Matrix.
   * Used by the orchestrator to synchronize damping values at boundaries after
   * initialization.
   */

  virtual VECTOR_REAL_VIEW& getDampingMatrix(int c) = 0;

  /**
   * @brief Access a Force Vector component.
   * Used by the orchestrator to synchronize forces at boundaries during each
   * time step.
   * @param component Component index (0 to getNumComponents()-1).
   */
  virtual VECTOR_REAL_VIEW& getForceVector(int component) = 0;

  /**
   * @brief Phase 1 of time step: Compute local forces (Stiffness + Source).
   *
   * This method computes the local contributions to the force vectors.
   * In a distributed simulation, the caller MUST synchronize the force vectors
   * via BoundarySynchronizer after calling this method and before calling
   * updateSolution.
   */
  virtual void computeForces(const float& dt, const int& timeSample,
                             DataStruct& data) = 0;

  /**
   * @brief Phase 2 of time step: Update solution using mass matrix and forces.
   *
   * Assumes force vectors have been fully assembled (including boundary
   * synchronization).
   */
  virtual void updateSolution(const float& dt, DataStruct& data) = 0;

  virtual void setAnisotropyType(model::AnisotropyType type) = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_SOLVER_H_
