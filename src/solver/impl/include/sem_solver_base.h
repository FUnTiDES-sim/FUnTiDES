//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.hpp: simple 2D acoustic wave equation solver
//
//  The SEMsolver class serves as a base class for the Spectral Element Method
//  solver. It provides core functionality to initialize FE operators,
//  advance pressure fields, apply forcing terms, and handle absorbing
//  boundaries.
//************************************************************************

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
   */
  virtual void computeFEInit(model::ModelApi<float, int> &mesh,
                             const std::array<float, 3> &sponge_size,
                             const bool surface_sponge,
                             const float taper_delta_) =0;


  /**
   * @brief Initialize arrays required by the finite element solver.
   */
  virtual void initFEarrays()=0;

  /**
   * @brief Allocate memory for FE-related arrays (mass, stiffness, etc.).
   */
  virtual void allocateFEarrays()=0;

  /**
   * @brief Initialize sponge (absorbing layer) coefficients.
   */
  virtual void initSpongeValues()=0;

  /**
   * @brief Reset global FE vectors (mass, stiffness) before accumulation.
   *
   * @param numNodes Total number of global nodes.
   */
  virtual void resetGlobalVectors(int numNodes)=0;

  /**
   * @brief Compute the global mass matrix.
   *  once at the beginning of the simulation.
   * @param isModelOnNodes True if the velocity model is defined on nodes, false
   *                       if on elements
   */
  virtual void computeGlobalMassMatrix()=0;

  virtual void outputSolutionValues(const int &indexTimeStep, int &i1,
                                    int &myElementSource,
                                    const ARRAY_REAL_VIEW &field, const char* fieldName)=0;

};

#endif  // SEM_SOLVERBASE_HPP_
