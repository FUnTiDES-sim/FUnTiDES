#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "data_type.h"

using namespace std::chrono;

struct SolverUtils
{
float evaluateRicker(float const& time_n, float const& f0, int order)
{
  float const o_tpeak = 1.0 / f0;  // tdelay dans Gar6more
  
  // Fenêtre temporelle : active seulement entre [0, 2*tdelay]
  if ((time_n <= 0.0f) || (time_n >= 2.0f * o_tpeak))
  {
    return 0.0f;
  }
  
  constexpr float pi = M_PI;
  float const alpha = -(f0 * pi) * (f0 * pi);  // alpha = -(π*f0)²
  float const time_d = time_n - o_tpeak;        // t1 dans Gar6more
  float const gaussian = exp(alpha * time_d * time_d);
  int const sgn = (order % 2 == 0) ? 1 : -1;    // (-1)^(deriv+1)
  
  float pulse = 0.0f;
  
  switch(order)
  {
    case 0:
      pulse = sgn * gaussian;
      break;
      
    case 1:
      pulse = sgn * (2.0f * alpha * time_d) * gaussian;
      break;
      
    case 2:
    {
      float const alpha2 = alpha * alpha;
      float const time_d2 = time_d * time_d;
      pulse = sgn * (2.0f * alpha + 4.0f * alpha2 * time_d2) * gaussian;
      break;
    }
    
    case 3:
    {
      float const alpha2 = alpha * alpha;
      float const alpha3 = alpha2 * alpha;
      float const time_d3 = time_d * time_d * time_d;
      pulse = sgn * (12.0f * alpha2 * time_d + 8.0f * alpha3 * time_d3) * gaussian;
      break;
    }
    
    case 4:
    {
      float const alpha2 = alpha * alpha;
      float const alpha3 = alpha2 * alpha;
      float const alpha4 = alpha3 * alpha;
      float const time_d2 = time_d * time_d;
      float const time_d4 = time_d2 * time_d2;
      pulse = sgn * (12.0f * alpha2 + 48.0f * alpha3 * time_d2 + 16.0f * alpha4 * time_d4) * gaussian;
      break;
    }
    
    default:
      std::cout << "rickerOrder must be between 0 and 4" << std::endl;
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
