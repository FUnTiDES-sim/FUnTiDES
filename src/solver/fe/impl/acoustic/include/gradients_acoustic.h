#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_GRADIENTS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_GRADIENTS_ACOUSTIC_H_
#include <data_type.h>

#include "gradients.h"

namespace solver
{
namespace fe
{
/**
 * @brief Acoustic gradients data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct GradientsAcoustic : public Gradients
{
  /// Number of gradients for model inversions (2 for Kappa and Buoyancy)
  static constexpr int kNumGrads = 2;

  /// Primary field name
  static constexpr const char* kGradsNames[2] = {"gradKappa","gradBuoyancy"};

  GradientsAcoustic(VECTOR_REAL_VIEW gradKappa,
                    VECTOR_REAL_VIEW gradBuoyancy)
      : m_gradKappa(gradKappa), m_gradBuoyancy(gradBuoyancy)
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
        return m_gradKappa;
      case 1:
        return m_gradBuoyancy;
      default:
        return m_gradKappa;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "Grad Kappa size: " << m_gradKappa.extent(0)
              << std::endl;
    std::cout << "Grad Buoyancy size: " << m_gradBuoyancy.extent(0)
              << std::endl;
  }

  VECTOR_REAL_VIEW m_gradKappa;  ///< Pressure field at previous time step
  VECTOR_REAL_VIEW m_gradBuoyancy;  ///< Pressure field at current time step
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_GRADIENTS_ACOUSTIC_H_
