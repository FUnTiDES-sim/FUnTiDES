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

#include "SolverBase.hpp"

#include "dataType.hpp"

#ifdef USE_KOKKOS
#include <KokkosExp_InterOp.hpp>
#endif

#include <cmath>
#include <model.hpp>


  struct SEMsolverData : SolverBase::DataStruct
  {
    SEMsolverData( int i1,
                   int i2,
                   ARRAY_REAL_VIEW const & rhsTerm,
                   ARRAY_REAL_VIEW const & pnGlobal,
                   VECTOR_INT_VIEW const & rhsElement,
                   ARRAY_REAL_VIEW const & rhsWeights):
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
          typename INTEGRAL_TYPE >
class SEMsolver : public SolverBase
{
public:

  constexpr static int numPoints1d = ORDER + 1;
  constexpr static int numPointsPerElem = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  SEMsolver() = default;
  ~SEMsolver() = default;



  /**
   * @brief Construct and initialize solver from a mesh.
   * @param mesh Reference to the mesh used in the simulation.
   */
  virtual void computeFEInit( Mesh const & mesh ) override final;

  /**
   * @brief Initialize all finite element structures:
   * basis functions, integrals, global arrays, etc.
   *
   * This function advances the pressure field `pnGlobal` by one time step using
   * a second-order explicit scheme. It resets global accumulators, applies RHS
   * forcing, computes local element contributions to mass and stiffness
   * matrices, updates pressure for interior nodes, and applies sponge damping.
   *
   * @param timeSample   Index of the current time step in `rhsTerm`
   * @param nPointsPerElement Number of quadrature points per element
   * @param i1           Index for pressure at previous time step
   * @param i2           Index for pressure at current time step
   * @param myInfo       Structure containing mesh and solver configuration
   * @param rhsTerm      External forcing term, function of space and time
   * @param pnGlobal     2D array storing the global pressure field [node][time]
   * @param rhsElement   List of elements with a non-zero forcing term
   */
  virtual void computeOneStep( const float & dt,
                              const int & timeSample,
                              DataStruct & data ) override final;

  virtual void outputPnValues( const int &indexTimeStep,
                               int &i1,
                               int &myElementSource,
                               const ARRAY_REAL_VIEW &pnGlobal ) override final;

  void initFEarrays();

  /**
   * @brief Initialize arrays required by the finite element solver.
   */
  void initSpongeValues(Mesh &mesh );

  // void spongeUpdate(const ARRAY_REAL_VIEW &pnGlobal, const int i1,
  //                   const int i2);

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
  void resetGlobalVectors();

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
  void applyRHSTerm(float const dt, int timeSample, int i2, const ARRAY_REAL_VIEW &rhsTerm,
                    const VECTOR_INT_VIEW &rhsElement,
                    const ARRAY_REAL_VIEW &pnGlobal,
                    const ARRAY_REAL_VIEW &rhsWeights);

  /**
   * @brief Assemble local element contributions to global FE vectors.
   *
   * @param i2       Current pressure field index
   * @param pnGlobal Global pressure field
   */
  void computeElementContributions( int nPointsPerElement,
                                    int i2,
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
  void updatePressureField( float const dt,
                     const ARRAY_REAL_VIEW &rhsTerm,
                     int i1,
                     int i2,
                     const ARRAY_REAL_VIEW &pnGlobal);


private:
  Mesh myMesh;                         ///< Internal copy of the mesh
  // const float myTimeStep = 0.001f;     ///< Time step size (fixed)
  float m_spongeSize = 250.;
  bool isSurface = true;

  Mesh m_mesh;
  // Basis functions and integrals
  INTEGRAL_TYPE myQkIntegrals;
  typename INTEGRAL_TYPE::PrecomputedData m_precomputedIntegralData;


  // Sponge tapering
  VECTOR_REAL_VIEW spongeTaperCoeff;

  VECTOR_REAL_VIEW quadraturePoints;
  VECTOR_REAL_VIEW weights;
  ARRAY_REAL_VIEW derivativeBasisFunction1D;

  // shared arrays
  VECTOR_REAL_VIEW massMatrixGlobal;
  VECTOR_REAL_VIEW yGlobal;
};

#endif // SEM_SOLVER_HPP_
