//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.hpp: simple 2D acoustic wave equation solver
//
//  The SEMsolver class serves as a base class for the Spectral Element Method
//  solver. It provides core functionality to initialize FE operators,
//  advance pressure fields, apply forcing terms, and handle absorbing boundaries.
//************************************************************************

#ifndef SEM_SOLVER_HPP_
#define SEM_SOLVER_HPP_

#include "dataType.hpp"
#include <BasisFunctions.hpp>
#include <Integrals.hpp>
#include <cmath>
#include <model.hpp>

#ifdef USE_KOKKOS
#include <KokkosExp_InterOp.hpp>
#endif

class SEMsolver {
public:
  /**
   * @brief Default constructor.
   */
  PROXY_HOST_DEVICE SEMsolver() {};

  /**
   * @brief Destructor.
   */
  PROXY_HOST_DEVICE ~SEMsolver(){};

  /**
   * @brief Construct and initialize solver from a mesh.
   * @param mesh Reference to the mesh used in the simulation.
   */
  SEMsolver(const Mesh mesh) { computeFEInit(mesh); };

  /**
   * @brief Initialize all finite element structures:
   * basis functions, integrals, global arrays, etc.
   *
   * @param mesh Mesh structure containing the domain information.
   */
  void computeFEInit(Mesh mesh);

  /**
   * @brief Compute one time step of the SEM wave equation solver.
   *
   * Advances the pressure field using explicit time integration.
   *
   * @param timeSample   Current time index into the RHS (source) term
   * @param dt           Delta time for this iteration
   * @param i1           Index for previous pressure field
   * @param i2           Index for current pressure field
   * @param rhsTerm      Right-hand side forcing term [node][time]
   * @param pnGlobal     Global pressure field [node][time]
   * @param rhsElement   List of active source elements
   * @param rhsWeights   Forcing weights per source node
   */
  void computeOneStep(const int &timeSample, const float dt,
                      const int &i1, const int &i2,
                      const ARRAY_REAL_VIEW &rhsTerm,
                      const ARRAY_REAL_VIEW &pnGlobal,
                      const VECTOR_INT_VIEW &rhsElement,
                      const ARRAY_REAL_VIEW &rhsWeights);

  /**
   * @brief Output pressure values at a specific time step.
   *
   * Typically used for recording seismograms or snapshots.
   *
   * @param mesh             Mesh structure
   * @param indexTimeStep    Time index to output
   * @param i1               Index for pressure buffer
   * @param myElementSource  Element containing the receiver
   * @param pnGlobal         Global pressure field [node][time]
   */
  void outputPnValues(Mesh mesh, const int &indexTimeStep, int &i1,
                      int &myElementSource, const ARRAY_REAL_VIEW &pnGlobal);

  /**
   * @brief Initialize arrays required by the finite element solver.
   */
  void initFEarrays();

  /**
   * @brief Allocate memory for FE-related arrays (mass, stiffness, etc.).
   */
  void allocateFEarrays();

  /**
   * @brief Initialize sponge (absorbing layer) coefficients.
   */
  void initSpongeValues();

  /**
   * @brief Reset global FE vectors (mass, stiffness) before accumulation.
   *
   * @param numNodes Total number of global nodes.
   */
  void resetGlobalVectors(int numNodes);

  /**
   * @brief Apply external forcing to the global pressure field.
   *
   * @param timeSample   Current time sample index
   * @param i2           Current pressure index
   * @param rhsTerm      RHS forcing term array
   * @param rhsElement   Indices of source elements
   * @param pnGlobal     Global pressure field (modified in-place)
   * @param rhsWeights   Forcing weights per node
   */
  void applyRHSTerm(int timeSample, float dt, int i2,
                    const ARRAY_REAL_VIEW &rhsTerm,
                    const VECTOR_INT_VIEW &rhsElement,
                    const ARRAY_REAL_VIEW &pnGlobal,
                    const ARRAY_REAL_VIEW &rhsWeights);

  /**
   * @brief Assemble local element contributions to global FE vectors.
   *
   * @param i2       Current pressure field index
   * @param pnGlobal Global pressure field
   */
  void computeElementContributions(int i2, const ARRAY_REAL_VIEW &pnGlobal);

  /**
   * @brief Update the global pressure field at interior nodes.
   *
   * Applies the time integration scheme.
   *
   * @param i1       Previous time step index
   * @param i2       Current time step index
   * @param pnGlobal Pressure field array (updated in-place)
   */
  void updatePressureField(float dt, int i1, int i2, const ARRAY_REAL_VIEW &pnGlobal);

  /**
   * @brief Accessor for the sponge tapering coefficients.
   * @return Const reference to the sponge taper coefficient vector.
   */
  const vectorReal &getSpongeTaperCoeff() const { return spongeTaperCoeff; }

private:
  Mesh myMesh;                         ///< Internal copy of the mesh
  // const float myTimeStep = 0.001f;     ///< Time step size (fixed)
  float m_spongeSize = 100.;
  bool isSurface = true;

  // Basis functions and integral objects
  SEMQkGLIntegrals myQkIntegrals;

  // Sponge tapering
  VECTOR_REAL_VIEW spongeTaperCoeff;

#ifdef USE_SEMCLASSIC
  SEMQkGLBasisFunctions myQkBasis;
  VECTOR_REAL_VIEW quadraturePoints;
  VECTOR_REAL_VIEW weights;
  ARRAY_REAL_VIEW derivativeBasisFunction1D;
#endif

  // Global FE vectors
  VECTOR_REAL_VIEW massMatrixGlobal;
  VECTOR_REAL_VIEW yGlobal;
};

#endif // SEM_SOLVER_HPP_
