
#ifndef SEM_SOLVERBASE_HPP_
#define SEM_SOLVERBASE_HPP_

#include <solver_base.h>

#include <cmath>

class SEMSolverBase : public SolverBase
{
 public:
  /**
   * @brief Initialize all finite element structures:
   * basis functions, integrals, global arrays, and sponge boundaries.
   *
   * @param mesh BaseMesh structure containing the domain information.
   * @param sponge_size Thickness (in elements) of absorbing sponge layers
   *                    in each direction [x, y, z] to prevent reflections.
   * @param surface_sponge Enable sponge at free surface (typically false
   *                       for geophysics to preserve natural reflections).
   * @param taper_delta_ Attenuation parameter for sponge layers.
   */
  virtual void computeFEInit(model::ModelApi<float, int> &mesh,
                             const std::array<float, 3> &sponge_size,
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
   * @brief Outputs solution field values at a specific time step
   *
   * This pure virtual function is responsible for writing or displaying
   * solution values for a given field at a particular time step. It must be
   * implemented by derived classes to define the specific output behavior.
   *
   * @param[in] indexTimeStep The index of the current time step for which
   *                          solution values are being output
   * @param[in,out] i1 Index variable (corresponding to the time n of the
   * solution)
   * @param[in,out] myElementSource Index or identifier of the source element
   *                                being processed
   * @param[in] field Constant view of the array containing the field values
   *                  to be output
   * @param[in] fieldName Name/identifier of the field being output (as a
   *                      C-string)
   */

  virtual void outputSolutionValues(const int &indexTimeStep, int &i1,
                                    int &myElementSource,
                                    const ARRAY_REAL_VIEW &field,
                                    const char *fieldName) = 0;
};

#endif  // SEM_SOLVERBASE_HPP_
