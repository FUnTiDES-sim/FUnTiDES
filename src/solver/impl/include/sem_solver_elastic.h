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

#ifndef SEM_SOLVER_ELASTIC_HPP_
#define SEM_SOLVER_ELASTIC_HPP_

#include <data_type.h>
#include <model.h>
#include <sem_solver_base.h>

#include <cmath>

struct SEMsolverDataElastic : public SolverBase::DataStruct
{
  SEMsolverDataElastic(int i1, int i2, ARRAY_REAL_VIEW rhsTermx,
                ARRAY_REAL_VIEW rhsTermy, ARRAY_REAL_VIEW rhsTermz,
                ARRAY_REAL_VIEW uxnGlobal,  ARRAY_REAL_VIEW uynGlobal, 
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
    std::cout << "SEMsolverDataElastic: i1=" << m_i1 << ", i2=" << m_i2 << std::endl;
    std::cout << "RHSx Term size: " << m_rhsTermx.extent(0) << std::endl;
    std::cout << "RHSy Term size: " << m_rhsTermy.extent(0) << std::endl;
    std::cout << "RHSz Term size: " << m_rhsTermz.extent(0) << std::endl;
    std::cout << "Uxn Global size: " << m_uxnGlobal.extent(0) << std::endl;
    std::cout << "Uyn Global size: " << m_uynGlobal.extent(0) << std::endl;
    std::cout << "Uzn Global size: " << m_uznGlobal.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_rhsElement.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_rhsWeights.extent(0) << std::endl;
  }

  int m_i1;
  int m_i2;
  ARRAY_REAL_VIEW m_rhsTermx;
  ARRAY_REAL_VIEW m_rhsTermy;
  ARRAY_REAL_VIEW m_rhsTermz;
  ARRAY_REAL_VIEW m_uxnGlobal;
  ARRAY_REAL_VIEW m_uynGlobal;
  ARRAY_REAL_VIEW m_uznGlobal;
  VECTOR_INT_VIEW m_rhsElement;
  ARRAY_REAL_VIEW m_rhsWeights;
};

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
class SEMsolverElastic : public SEMsolverBase
{
 public:
  SEMsolverElastic() = default;
  
  ~SEMsolverElastic() = default;

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

  void outputUnValues(const int &indexTimeStep, int &i1,
                      int &myElementSource,
                      const ARRAY_REAL_VIEW &uxnGlobal,
                      const ARRAY_REAL_VIEW &uynGlobal,
                      const ARRAY_REAL_VIEW &uznGlobal);

  void applyRHSTerm(int timeSample, float dt, int i2,
                    const ARRAY_REAL_VIEW &rhsTermx,
                    const ARRAY_REAL_VIEW &rhsTermy,
                    const ARRAY_REAL_VIEW &rhsTermz,
                    const VECTOR_INT_VIEW &rhsElement,
                    const ARRAY_REAL_VIEW &pnGlobal,
                    const ARRAY_REAL_VIEW &rhsWeights);

  void computeElementContributions(int i2, const ARRAY_REAL_VIEW &pnGlobal,
                                   bool isModelOnNodes);

  void updateDisplacementField(float dt, int i1, int i2,
                               const ARRAY_REAL_VIEW &uxnGlobal,
                               const ARRAY_REAL_VIEW &uynGlobal,
                               const ARRAY_REAL_VIEW &uznGlobal);

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
  VECTOR_REAL_VIEW uxGlobal;
  VECTOR_REAL_VIEW uyGlobal;
  VECTOR_REAL_VIEW uzGlobal;
};

#endif  // SEM_SOLVER_ELASTIC_HPP_
