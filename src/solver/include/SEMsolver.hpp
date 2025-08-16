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

#include "cartesian_sem_mesh.h"
#include "dataType.hpp"
#include "SolverBase.hpp"
#include <cmath>
#include <model.hpp>


struct SEMsolverData : SolverBase::DataStruct
{
  SEMsolverData( int i1,
                 int i2,
                 ARRAY_REAL_VIEW const & rhsTerm,
                 ARRAY_REAL_VIEW const & pnGlobal,
                 VECTOR_INT_VIEW const & rhsElement,
                 ARRAY_REAL_VIEW & rhsWeights):
    m_i1(i1),
    m_i2(i2),
    m_rhsTerm(rhsTerm),
    m_pnGlobal(pnGlobal),
    m_rhsElement(rhsElement),
    m_rhsWeights(rhsWeights)
  {}

  int m_i1;
  int m_i2;
  ARRAY_REAL_VIEW const & m_rhsTerm;
  ARRAY_REAL_VIEW const & m_pnGlobal;
  VECTOR_INT_VIEW const & m_rhsElement;
  ARRAY_REAL_VIEW const & m_rhsWeights;
};


template< int ORDER,
          typename INTEGRAL_TYPE,
          typename MESH_TYPE>
class SEMsolver : public SolverBase
{
public:
  /**
   * @brief Default constructor.
   */
  SEMsolver() = default;

  /**
   * @brief Destructor.
   */
  ~SEMsolver() = default;

  /**
   * @brief Initialize all finite element structures:
   * basis functions, integrals, global arrays, etc.
   *
   * @param mesh BaseMesh structure containing the domain information.
   */
  virtual void computeFEInit( BaseMesh<discretization_t, index_t> & mesh ) override final;

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
  virtual void computeOneStep( const float & dt,
                               const int & timeSample,
                               DataStruct & data ) override final;

  /**
   * @brief Output pressure values at a specific time step.
   *
   * Typically used for recording seismograms or snapshots.
   *
   * @param indexTimeStep    Time index to output
   * @param i1               Index for pressure buffer
   * @param myElementSource  Element containing the receiver
   * @param pnGlobal         Global pressure field [node][time]
   */
  virtual void outputPnValues( const int &indexTimeStep,
                               int &i1,
                               int &myElementSource,
                               const ARRAY_REAL_VIEW &pnGlobal) override final;

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
  void computeElementContributions(int i2,
                                   const ARRAY_REAL_VIEW &pnGlobal);

  /**
   * @brief Update the global pressure field at interior nodes.
   *
   * Applies the time integration scheme.
   *
   * @param i1       Previous time step index
   * @param i2       Current time step index
   * @param pnGlobal Pressure field array (updated in-place)
   */
  void updatePressureField(float dt,
                           int i1,
                           int i2,
                           const ARRAY_REAL_VIEW &pnGlobal);


private:
  MESH_TYPE m_mesh;

  static constexpr int nPointsElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float m_spongeSize = 250.;
  bool isSurface = true;

  // Basis functions and integral objects
  INTEGRAL_TYPE myQkIntegrals;
  typename INTEGRAL_TYPE::PrecomputedData m_precomputedIntegralData;

  // Sponge tapering
  VECTOR_REAL_VIEW spongeTaperCoeff;

  // Global FE vectors
  VECTOR_REAL_VIEW massMatrixGlobal;
  VECTOR_REAL_VIEW yGlobal;

  void computeFEInit(CartesianSEMmesh<discretization_t, index_t, ORDER> const & mesh);
};

#endif // SEM_SOLVER_HPP_
