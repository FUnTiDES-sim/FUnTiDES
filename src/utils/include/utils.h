#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "data_type.h"

using namespace std::chrono;

struct SolverUtils
{
  float evaluateRicker(float const& time_n, float const& f0, int order)
  {
    float const o_tpeak = 1.0 / f0;
    float pulse = 0.0;
    if ((time_n <= -0.9 * o_tpeak) || (time_n >= 2.9 * o_tpeak))
    {
      return pulse;
    }

    constexpr float pi = M_PI;
    float const lam = (f0 * pi) * (f0 * pi);

    switch (order)
    {
      case 2: {
        pulse = 2.0 * lam *
                (2.0 * lam * (time_n - o_tpeak) * (time_n - o_tpeak) - 1.0) *
                exp(-lam * (time_n - o_tpeak) * (time_n - o_tpeak));
      }
      break;
      case 1: {
        pulse = -2.0 * lam * (time_n - o_tpeak) *
                exp(-lam * (time_n - o_tpeak) * (time_n - o_tpeak));
      }
      break;
      case 0: {
        pulse = -(time_n - o_tpeak) *
                exp(-2 * lam * (time_n - o_tpeak) * (time_n - o_tpeak));
      }
      break;
      default:
        std::cout
            << "This option is not supported yet, rickerOrder must be 0, 1 or 2"
            << std::endl;
        break;
    }

    return pulse;
  }

  std::vector<float> computeSourceTerm(const int nSamples,
                                       const float timeSample, const float f0,
                                       const int order)
  {
    std::vector<float> sourceTerm(nSamples);
    for (int i = 0; i < nSamples; i++)
    {
      float time_n = i * timeSample;
      sourceTerm[i] = evaluateRicker(time_n, f0, order);
    }
    return sourceTerm;
  }
};
#endif  // UTILS_HPP_
