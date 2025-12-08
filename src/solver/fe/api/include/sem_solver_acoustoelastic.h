#ifndef SEM_SOLVER_ACOUSTOELASTIC_HPP_
#define SEM_SOLVER_ACOUSTOELASTIC_HPP_

#include <sem_solver_acoustic.h>
#include <sem_solver_elastic.h>

// Element type constants
constexpr int ELEMENT_TYPE_ACOUSTIC = 1;
constexpr int ELEMENT_TYPE_ELASTIC = 2;

/**
 * @brief Data structure for acoustoelastic wave propagation solver.
 *
 * Uses 3 time buffers for elastic fields (n-1, n, n+1) to compute acceleration for coupling.
 * Acoustic fields use standard 2 buffers.
 */
struct SEMsolverDataAcoustoElastic : public SolverBase::DataStruct
{
  SEMsolverDataAcoustoElastic(int i0, int i1, int i2,
                              ARRAY_REAL_VIEW rhsTerm,
                              ARRAY_REAL_VIEW rhsTermx,
                              ARRAY_REAL_VIEW rhsTermy,
                              ARRAY_REAL_VIEW rhsTermz,
                              ARRAY_REAL_VIEW pnGlobal,
                              ARRAY_REAL_VIEW uxnGlobal,
                              ARRAY_REAL_VIEW uynGlobal,
                              ARRAY_REAL_VIEW uznGlobal,
                              VECTOR_INT_VIEW rhsElement,
                              ARRAY_REAL_VIEW rhsWeights)
      : m_i0(i0),
        m_i1(i1),
        m_i2(i2),
        m_rhsTerm(rhsTerm),
        m_rhsTermx(rhsTermx),
        m_rhsTermy(rhsTermy),
        m_rhsTermz(rhsTermz),
        m_pnGlobal(pnGlobal),
        m_uxnGlobal(uxnGlobal),
        m_uynGlobal(uynGlobal),
        m_uznGlobal(uznGlobal),
        m_rhsElement(rhsElement),
        m_rhsWeights(rhsWeights)
  {
  }

  void print() const override
  {
    std::cout << "SEMsolverDataAcoustoElastic: i0=" << m_i0 << ", i1=" << m_i1 << ", i2=" << m_i2 << std::endl;
    std::cout << "Acoustic RHS size: " << m_rhsTerm.extent(0) << std::endl;
    std::cout << "Elastic RHSx size: " << m_rhsTermx.extent(0) << std::endl;
    std::cout << "Pressure field size: " << m_pnGlobal.extent(0) << " x " << m_pnGlobal.extent(1) << std::endl;
    std::cout << "Displacement Ux size: " << m_uxnGlobal.extent(0) << " x " << m_uxnGlobal.extent(1) << std::endl;
  }

  int m_i0;  ///< Oldest time index (n-1 for elastic)
  int m_i1;  ///< Middle time index (n for elastic, n-1 for acoustic)
  int m_i2;  ///< Newest time index (n+1 for elastic, n for acoustic)
  
  ARRAY_REAL_VIEW m_rhsTerm;     ///< Acoustic RHS forcing term
  ARRAY_REAL_VIEW m_rhsTermx;    ///< Elastic X-component forcing term
  ARRAY_REAL_VIEW m_rhsTermy;    ///< Elastic Y-component forcing term
  ARRAY_REAL_VIEW m_rhsTermz;    ///< Elastic Z-component forcing term
  ARRAY_REAL_VIEW m_pnGlobal;    ///< Acoustic pressure field (numNodes, 2)
  ARRAY_REAL_VIEW m_uxnGlobal;   ///< Elastic X-displacement field (numNodes, 3)
  ARRAY_REAL_VIEW m_uynGlobal;   ///< Elastic Y-displacement field (numNodes, 3)
  ARRAY_REAL_VIEW m_uznGlobal;   ///< Elastic Z-displacement field (numNodes, 3)
  VECTOR_INT_VIEW m_rhsElement;  ///< Source element indices
  ARRAY_REAL_VIEW m_rhsWeights;  ///< Forcing weights per node
};

/**
 * @brief Spectral Element Method solver for acoustoelastic wave propagation.
 * 
 * Uses composition to combine acoustic and elastic solvers, avoiding diamond inheritance.
 * Handles coupling at the interface between acoustic (fluid) and elastic (solid) media.
 * Elements are tagged with integer types: 1=acoustic, 2=elastic.
 * 
 * Uses 3 time buffers for elastic displacement to compute acceleration for coupling.
 * 
 * SIMPLIFIED VERSION: Assumes horizontal interface with normal = (0, 0, 1).
 * TODO: Implement proper interface geometry computation with face normals and areas.
 *
 * @tparam ORDER Polynomial order of the spectral elements
 * @tparam INTEGRAL_TYPE Type for numerical integration (basis functions, quadrature)
 * @tparam MESH_TYPE Type of the computational mesh
 * @tparam IS_MODEL_ON_NODES Boolean to say if the model is located on nodes (true) or on elements (false)
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class SEMsolverAcoustoElastic : public SEMSolverBase
{
 public:
  SEMsolverAcoustoElastic() = default;
  ~SEMsolverAcoustoElastic() = default;

  void computeFEInit(model::ModelApi<float, int> &mesh,
                     const std::array<float, 3> &sponge_size,
                     const bool surface_sponge,
                     const float taper_delta_) override;

  void computeOneStep(const float &dt, const int &timeSample,
                      DataStruct &data) override;

  void initFEarrays() override;
  void allocateFEarrays() override;
  void initSpongeValues() override;
  void resetGlobalVectors(int numNodes) override;
  void computeGlobalMassMatrix() override;

  void outputSolutionValues(const int &indexTimeStep, int &i1,
                            int &myElementSource, const ARRAY_REAL_VIEW &field,
                            const char *fieldName) override;

  void tagElements();
  void tagNodes();

  int getElementType(int elementIdx) const { return m_elementType(elementIdx); }
  int getNumAcousticElements() const { return m_numAcousticElements; }
  int getNumElasticElements() const { return m_numElasticElements; }
  int getNumInterfaceNodes() const { return m_numInterfaceNodes; }

 private:
  /**
   * @brief Apply acoustic-to-elastic coupling.
   * 
   * Adds pressure force to elastic displacement: u_n+1 += dt² * (-p_n) * normal / M_elastic
   * Called after elastic solver completes.
   * 
   * SIMPLIFIED: Uses vertical normal (0, 0, 1) for horizontal interface.
   */
  void applyCouplingAcousticToElastic(float dt, SEMsolverDataAcoustoElastic &data);

  /**
   * @brief Apply elastic-to-acoustic coupling.
   * 
   * Adds elastic acceleration to acoustic pressure: p_n+1 += accel_elastic * normal / M_acoustic
   * Called after acoustic solver completes.
   * Uses 3 time levels (i0, i1, i2) to compute acceleration.
   * 
   * SIMPLIFIED: Uses vertical normal (0, 0, 1) for horizontal interface.
   */
  void applyCouplingElasticToAcoustic(float dt, SEMsolverDataAcoustoElastic &data);

  void identifyInterfaceElements();

 private:
  SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES> 
      m_acousticSolver;
  SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES> 
      m_elasticSolver;

  MESH_TYPE m_mesh;

  VECTOR_INT_VIEW m_elementType;        ///< Element type (1=acoustic, 2=elastic)
  VECTOR_BOOL_VIEW m_isInterfaceNode;   ///< True if node is at interface

  int m_numAcousticElements{0};
  int m_numElasticElements{0};
  int m_numInterfaceNodes{0};

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;
  
  static constexpr float mu_tolerance_ = 1.0e-6f;
};

#endif  // SEM_SOLVER_ACOUSTOELASTIC_HPP_
