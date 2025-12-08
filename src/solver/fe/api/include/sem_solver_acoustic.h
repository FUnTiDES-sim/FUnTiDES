#ifndef SEM_SOLVER_ACOUSTIC_HPP_
#define SEM_SOLVER_ACOUSTIC_HPP_

#include <data_type.h>
#include <model.h>
#include <sem_solver_base.h>

#include <cmath>

struct SEMsolverDataAcoustic : public SolverBase::DataStruct
{
  SEMsolverDataAcoustic(int i1, int i2, ARRAY_REAL_VIEW rhsTerm,
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

  int m_i1;                      ///< Previous time step index
  int m_i2;                      ///< Current time step index
  ARRAY_REAL_VIEW m_rhsTerm;     ///< RHS forcing term
  ARRAY_REAL_VIEW m_pnGlobal;    ///< Pressure field
  VECTOR_INT_VIEW m_rhsElement;  ///< Source element indices
  ARRAY_REAL_VIEW m_rhsWeights;  ///< Forcing weights per node
};

/**
 * @brief Spectral Element Method solver for acoustic wave propagation.
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
class SEMsolverAcoustic : public SEMSolverBase
{
 public:
  SEMsolverAcoustic() = default;
  ~SEMsolverAcoustic() = default;

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
   * @brief Compute one time step of the acoustic wave equation solver.
   *
   * Advances the pressure field using explicit time integration.
   * If an element type mask is set, only processes elements of the specified type.
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
   * @param fieldName The name of the field to output (here it can be pnGlobal)
   */
  void outputSolutionValues(const int &indexTimeStep, int &i1,
                            int &myElementSource, const ARRAY_REAL_VIEW &field,
                            const char *fieldName) override;

  /**
   * @brief Apply external forcing to the global pressure field.
   *
   * @param timeSample Current time sample index.
   * @param dt Delta time for this iteration.
   * @param i2 Current displacement index.
   * @param rhsTerm RHS forcing term array.
   * @param rhsElement Indices of source elements.
   * @param rhsWeights Forcing weights per node.
   */
  void applyRHSTerm(int timeSample, float dt, int i2,
                    const ARRAY_REAL_VIEW &rhsTerm,
                    const VECTOR_INT_VIEW &rhsElement,
                    const ARRAY_REAL_VIEW &rhsWeights);

  /**
   * @brief Assemble local element contributions to global FE vectors.
   *
   * @param i2 Current displacement field index.
   * @param pnGlobal Global pressure field.
   */
  void computeElementContributions(int i2, const ARRAY_REAL_VIEW &pnGlobal);

  /**
   * @brief Update the global pressure field at interior nodes.
   *
   * Applies the time integration scheme for acoustic wave propagation.
   *
   * @param dt Delta time for this iteration.
   * @param i1 Previous time step index.
   * @param i2 Current time step index.
   * @param pnGlobal Pressure field array (updated in-place).
   */
  void updatePressureField(float dt, int i1, int i2,
                           const ARRAY_REAL_VIEW &pnGlobal);

  /**
   * @brief Set element type mask for selective computation.
   *
   * When set, only elements where elementType[i] == myType are processed.
   * If nullptr (default), all elements are processed.
   * Used by acoustoelastic solver to restrict acoustic solver to fluid regions.
   *
   * @param elementType Pointer to integer array (one entry per element, values: 1=acoustic, 2=elastic).
   * @param myType Element type to process (1 for acoustic).
   */
  void setElementMask(const VECTOR_INT_VIEW *elementType, int myType) {
    m_elementType = elementType;
    m_myType = myType;
  }

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
  VECTOR_REAL_VIEW yGlobal;           ///< Global pressure work vector

  /// Element type mask for selective computation (nullptr = process all elements)
  const VECTOR_INT_VIEW *m_elementType{nullptr};
  int m_myType{0};  ///< Element type this solver should process (1 for acoustic)
};

#endif  // SEM_SOLVER_ACOUSTIC_HPP_
