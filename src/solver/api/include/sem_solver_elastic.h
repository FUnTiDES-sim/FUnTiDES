#ifndef SEM_SOLVER_ELASTIC_HPP_
#define SEM_SOLVER_ELASTIC_HPP_

#include <data_type.h>
#include <model.h>
#include <sem_solver_base.h>

#include <cmath>

/**
 * @brief Data structure for elastic wave propagation solver.
 *
 * Contains displacement fields, forcing terms, and time indices
 * for elastic wave equation computations.
 */
struct SEMsolverDataElastic : public SolverBase::DataStruct
{
  SEMsolverDataElastic(int i1, int i2, ARRAY_REAL_VIEW rhsTermx,
                       ARRAY_REAL_VIEW rhsTermy, ARRAY_REAL_VIEW rhsTermz,
                       ARRAY_REAL_VIEW uxnGlobal, ARRAY_REAL_VIEW uynGlobal,
                       ARRAY_REAL_VIEW uznGlobal, VECTOR_INT_VIEW rhsElement,
                       ARRAY_REAL_VIEW rhsWeights)
      : m_i1(i1),
        m_i2(i2),
        m_rhsTermx(rhsTermx),
        m_rhsTermy(rhsTermy),
        m_rhsTermz(rhsTermz),
        m_uxnGlobal(uxnGlobal),
        m_uynGlobal(uynGlobal),
        m_uznGlobal(uznGlobal),
        m_rhsElement(rhsElement),
        m_rhsWeights(rhsWeights)
  {
  }

  void print() const override
  {
    std::cout << "SEMsolverDataElastic: i1=" << m_i1 << ", i2=" << m_i2
              << std::endl;
    std::cout << "RHSx Term size: " << m_rhsTermx.extent(0) << std::endl;
    std::cout << "RHSy Term size: " << m_rhsTermy.extent(0) << std::endl;
    std::cout << "RHSz Term size: " << m_rhsTermz.extent(0) << std::endl;
    std::cout << "Uxn Global size: " << m_uxnGlobal.extent(0) << std::endl;
    std::cout << "Uyn Global size: " << m_uynGlobal.extent(0) << std::endl;
    std::cout << "Uzn Global size: " << m_uznGlobal.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_rhsElement.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_rhsWeights.extent(0) << std::endl;
  }

  int m_i1;                      ///< Previous time step index
  int m_i2;                      ///< Current time step index
  ARRAY_REAL_VIEW m_rhsTermx;    ///< X-component forcing term
  ARRAY_REAL_VIEW m_rhsTermy;    ///< Y-component forcing term
  ARRAY_REAL_VIEW m_rhsTermz;    ///< Z-component forcing term
  ARRAY_REAL_VIEW m_uxnGlobal;   ///< X-displacement field
  ARRAY_REAL_VIEW m_uynGlobal;   ///< Y-displacement field
  ARRAY_REAL_VIEW m_uznGlobal;   ///< Z-displacement field
  VECTOR_INT_VIEW m_rhsElement;  ///< Source element indices
  ARRAY_REAL_VIEW m_rhsWeights;  ///< Forcing weights per node
};

/**
 * @brief Spectral Element Method solver for elastic wave propagation.
 *
 * @tparam ORDER Polynomial order of the spectral elements
 * @tparam INTEGRAL_TYPE Type for numerical integration (basis functions,
 * quadrature)
 * @tparam MESH_TYPE Type of the computational mesh
 * @tparam IS_MODEL_ON_NODES Boolean to say if the model is located on nodes
 * (true) or on elements (false)
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class SEMsolverElastic : public SEMSolverBase
{
 public:
  /**
   * @brief Default constructor.
   */
  SEMsolverElastic() = default;

  /**
   * @brief Destructor.
   */
  ~SEMsolverElastic() = default;

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
  void computeFEInit(model::ModelApi<float, int> &mesh,
                     const std::array<float, 3> &sponge_size,
                     const bool surface_sponge,
                     const float taper_delta_) override;

  /**
   * @brief Compute one time step of the elastic wave equation solver.
   *
   * Advances the displacement field using explicit time integration.
   *
   * @param dt Delta time for this iteration.
   * @param timeSample Current time index into the RHS (source) term.
   * @param data DataStruct containing all necessary arrays.
   */
  void computeOneStep(const float &dt, const int &timeSample,
                      DataStruct &data) override;

  /**
   * @brief Initialize arrays required by the finite element solver.
   */
  void initFEarrays() override;

  /**
   * @brief Allocate memory for FE-related arrays (mass, stiffness, etc.).
   */
  void allocateFEarrays() override;

  /**
   * @brief Initialize sponge (absorbing layer) coefficients.
   */
  void initSpongeValues() override;

  /**
   * @brief Reset global FE vectors (mass, stiffness) before accumulation.
   *
   * @param numNodes Total number of global nodes.
   */
  void resetGlobalVectors(int numNodes) override;

  /**
   * @brief Compute the global mass matrix, accounting for the model.
   */
  void computeGlobalMassMatrix() override;

  /**
   * @brief Output displacement values at a specific time step.
   *
   * Typically used for recording seismograms or snapshots.
   *
   * @param indexTimeStep Time index to output.
   * @param i1 Index for displacement buffer.
   * @param myElementSource Element containing the receiver.
   * @param field The field to output
   * @param fieldName The name of the field to output (here it can be uxnGlobal,
   * uynGloba, uznGlobal)
   */
  void outputSolutionValues(const int &indexTimeStep, int &i1,
                            int &myElementSource,
                            const ARRAY_REAL_VIEW &uxnGlobal,
                            const char *fieldName) override;

  /**
   * @brief Apply external forcing to the global displacement field.
   *
   * @param timeSample Current time sample index.
   * @param dt Delta time for this iteration.
   * @param i2 Current displacement index.
   * @param rhsTermx X-component RHS forcing term array.
   * @param rhsTermy Y-component RHS forcing term array.
   * @param rhsTermz Z-component RHS forcing term array.
   * @param rhsElement Indices of source elements.
   * @param rhsWeights Forcing weights per node.
   */
  void applyRHSTerm(int timeSample, float dt, int i2,
                    const ARRAY_REAL_VIEW &rhsTermx,
                    const ARRAY_REAL_VIEW &rhsTermy,
                    const ARRAY_REAL_VIEW &rhsTermz,
                    const VECTOR_INT_VIEW &rhsElement,
                    const ARRAY_REAL_VIEW &rhsWeights);

  /**
   * @brief Assemble local element contributions to global FE vectors.
   *
   * @param i2 Current displacement field index.
   * @param uxnGlobal Global displacement in x direction  field.
   * @param uynGlobal Global displacement in y direction  field.
   * @param uznGlobal Global displacement in z direction  field.
   */
  void computeElementContributions(int i2, const ARRAY_REAL_VIEW &uxnGlobal,
                                   const ARRAY_REAL_VIEW &uynGlobal,
                                   const ARRAY_REAL_VIEW &uznGlobal);

  /**
   * @brief Update the global displacement field at interior nodes.
   *
   * Applies the time integration scheme for elastic wave propagation.
   *
   * @param dt Delta time for this iteration.
   * @param i1 Previous time step index.
   * @param i2 Current time step index.
   * @param uxnGlobal X-displacement field array (updated in-place).
   * @param uynGlobal Y-displacement field array (updated in-place).
   * @param uznGlobal Z-displacement field array (updated in-place).
   */
  void updateDisplacementField(float dt, int i1, int i2,
                               const ARRAY_REAL_VIEW &uxnGlobal,
                               const ARRAY_REAL_VIEW &uynGlobal,
                               const ARRAY_REAL_VIEW &uznGlobal);

  /**
   * @brief Compute the elasticity matrix at a given node.
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
  PROXY_HOST_DEVICE
  void computeCMatrix(float const vp, float const vs, float const rho,
                      float const delta, float const epsilon, float const gamma,
                      float const phi, float const theta,
                      float (&C)[6][6]) const;

 private:
  MESH_TYPE m_mesh;  ///< Computational mesh

  /// Number of nodes per element
  static constexpr int nPointsElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float sponge_size_[3];  ///< Sponge layer thickness [x, y, z]
  bool surface_sponge_;   ///< Enable sponge at free surface
  float taper_delta_;     ///< Attenuation parameter

  /// Basis functions and integral objects
  INTEGRAL_TYPE myQkIntegrals;

  VECTOR_REAL_VIEW spongeTaperCoeff;  ///< Sponge tapering coefficients
  VECTOR_REAL_VIEW massMatrixGlobal;  ///< Global mass matrix
  VECTOR_REAL_VIEW uxGlobal;  ///< Global displacment X-component work vector
  VECTOR_REAL_VIEW uyGlobal;  ///< Global displacment Y-component work vector
  VECTOR_REAL_VIEW uzGlobal;  ///< Global displacment Z-component work vector
};

#endif  // SEM_SOLVER_ELASTIC_HPP_
