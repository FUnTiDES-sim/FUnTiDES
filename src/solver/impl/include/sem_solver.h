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

#ifndef SEM_SOLVER_HPP_
#define SEM_SOLVER_HPP_

#include <data_type.h>
#include <model.h>
#include <solver_base.h>

#include <cmath>

struct SEMsolverData : SolverBase::DataStruct
{
  SEMsolverData(int i1, int i2, ARRAY_REAL_VIEW rhsTerm,
                ARRAY_REAL_VIEW pnGlobal, VECTOR_INT_VIEW rhsElement,
                ARRAY_REAL_VIEW rhsWeights)
      : m_i1(i1),
        m_i2(i2),
        m_rhsTerm(rhsTerm),
        m_pnGlobal(pnGlobal),
        m_rhsElement(rhsElement),
        m_rhsWeights(rhsWeights)
  {
  }

  void print() const override
  {
    std::cout << "SEMsolverData: i1=" << m_i1 << ", i2=" << m_i2 << std::endl;
    std::cout << "RHS Term size: " << m_rhsTerm.extent(0) << std::endl;
    std::cout << "Pn Global size: " << m_pnGlobal.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_rhsElement.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_rhsWeights.extent(0) << std::endl;
  }

  int m_i1;
  int m_i2;
  ARRAY_REAL_VIEW m_rhsTerm;
  ARRAY_REAL_VIEW m_pnGlobal;
  VECTOR_INT_VIEW m_rhsElement;
  ARRAY_REAL_VIEW m_rhsWeights;
};

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
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
                             const float taper_delta_) override final;
  /**
   * @brief Compute one time step of the SEM wave equation solver.
   *
   * Advances the pressure field using explicit time integration.
   *
   * @param timeSample   Current time index into the RHS (source) term
   * @param dt           Delta time for this iteration
   * @param data         DataStruct containing all necessary arrays
   * @param isModelOnNodes True if the velocity model is defined on nodes, false
   * if on elements
   */
  virtual void computeOneStep(const float &dt, const int &timeSample,
                              DataStruct &data) override final;

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
  virtual void outputPnValues(const int &indexTimeStep, int &i1,
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
   * @brief Compute the global mass matrix, accounting for the model.
   *
   * @param isModelOnNodes True if the velocity model is defined on nodes, false
   *                       if on elements
   */
  void computeGlobalMassMatrix(bool isModelOnNodes);

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
  void computeElementContributions(int i2, const ARRAY_REAL_VIEW &pnGlobal,
                                   bool isModelOnNodes);

  /**
   * @brief Update the global pressure field at interior nodes.
   *
   * Applies the time integration scheme.
   *
   * @param i1       Previous time step index
   * @param i2       Current time step index
   * @param pnGlobal Pressure field array (updated in-place)
   */
  void updatePressureField(float dt, int i1, int i2,
                           const ARRAY_REAL_VIEW &pnGlobal);

 private:
  MESH_TYPE m_mesh;
  const MESH_TYPE *m_mesh_ptr;

  static constexpr int nPointsElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;  // attenuation parameter

  // Basis functions and integral objects
  INTEGRAL_TYPE myQkIntegrals;
  typename INTEGRAL_TYPE::PrecomputedData m_precomputedIntegralData;

  // Sponge tapering
  VECTOR_REAL_VIEW spongeTaperCoeff;

  // Global FE vectors
  VECTOR_REAL_VIEW massMatrixGlobal;
  VECTOR_REAL_VIEW yGlobal;

  void computeFEInit(MESH_TYPE const &mesh);
};

#endif  // SEM_SOLVER_HPP_
