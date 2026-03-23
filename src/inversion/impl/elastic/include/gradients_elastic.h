#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_GRADIENTS_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_GRADIENTS_ELASTIC_H_
#include <data_type.h>

#include "gradients.h"

namespace solver
{
namespace fe
{
/**
 * @brief Elastic gradients data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct GradientsElastic : public Gradients
{
  /// Number of gradients for model inversions (3 for Rho, Lambda and Mu)
  static constexpr int kNumGrads = 3;

  /// Primary field name
  static constexpr const char* kGradsNames[3] = {"gradRho","gradLambda","gradMu"};

  GradientsElastic(VECTOR_REAL_VIEW gradRho,
                   VECTOR_REAL_VIEW gradLambda,
                   VECTOR_REAL_VIEW gradMu)
      : m_gradRho(gradRho), m_gradLambda(gradLambda), m_gradMu(gradMu)
  {
  }

  int getNumGrads() const override final { return kNumGrads; }

  const char* const* getGradsNames() const override final
  {
    return kGradsNames;
  }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getCurrentGrads(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_gradRho;
      case 1:
        return m_gradLambda;
      case 2:
        return m_gradMu;
      default:
        return m_gradRho;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "Grad Rho size: " << m_gradRho.extent(0)
              << std::endl;
    std::cout << "Grad Lambda size: " << m_gradLambda.extent(0)
              << std::endl;
    std::cout << "Grad Mu size: " << m_gradMu.extent(0)
              << std::endl;
  }

  VECTOR_REAL_VIEW m_gradRho;     ///< Gradient Rho field
  VECTOR_REAL_VIEW m_gradLambda;  ///< Gradient Lambda field
  VECTOR_REAL_VIEW m_gradMu;      ///< Gradient Mu field
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_GRADIENTS_ELASTIC_H_
