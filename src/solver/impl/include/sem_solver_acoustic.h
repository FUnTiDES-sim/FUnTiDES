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

  int m_i1;
  int m_i2;
  ARRAY_REAL_VIEW m_rhsTerm;
  ARRAY_REAL_VIEW m_pnGlobal;
  VECTOR_INT_VIEW m_rhsElement;
  ARRAY_REAL_VIEW m_rhsWeights;
};

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class SEMsolverAcoustic : public SEMSolverBase
{
 public:
  SEMsolverAcoustic() = default;
  ~SEMsolverAcoustic() = default;

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
                            int &myElementSource,
                            const ARRAY_REAL_VIEW &pnGlobal,
                            const char *fieldName) override;

  void applyRHSTerm(int timeSample, float dt, int i2,
                    const ARRAY_REAL_VIEW &rhsTerm,
                    const VECTOR_INT_VIEW &rhsElement,
                    const ARRAY_REAL_VIEW &pnGlobal,
                    const ARRAY_REAL_VIEW &rhsWeights);

  void computeElementContributions(int i2, const ARRAY_REAL_VIEW &pnGlobal);

  void updatePressureField(float dt, int i1, int i2,
                           const ARRAY_REAL_VIEW &pnGlobal);

 private:
  MESH_TYPE m_mesh;

  static constexpr int nPointsElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float sponge_size_[3];
  bool surface_sponge_;
  float taper_delta_;

  INTEGRAL_TYPE myQkIntegrals;
  typename INTEGRAL_TYPE::PrecomputedData m_precomputedIntegralData;

  VECTOR_REAL_VIEW spongeTaperCoeff;
  VECTOR_REAL_VIEW massMatrixGlobal;
  VECTOR_REAL_VIEW yGlobal;
};

#endif  // SEM_SOLVER_ACOUSTIC_HPP_
