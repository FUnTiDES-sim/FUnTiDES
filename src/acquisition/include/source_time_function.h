#ifndef FUNTIDES_ACQUISITION_INCLUDE_SOURCE_TIME_FUNCTION_H_
#define FUNTIDES_ACQUISITION_INCLUDE_SOURCE_TIME_FUNCTION_H_
#include "data_type.h"

using namespace std::chrono;

/**
 * @brief Evaluates source time functions of the Ricker family and samples them on a regular time grid.
 *
 * The functions are stateless. Time and peak time are expected in the same unit as 1/f0.
 * @todo VERIFY: is the time unit seconds and f0 in Hz?
 */
struct SourceTimeFunction {
  /**
   * @brief Evaluates the source time function of a given order at one time.
   *
   * The value is forced to 0 outside the open window (-0.9 * tpeak, 2.9 * tpeak).
   * An unsupported order prints a message to std::cout and yields 0.
   *
   * @param[in] time_n Evaluation time.
   * @param[in] f0 Central frequency.
   * @param[in] order Function selector, 0 to 4.
   * @todo VERIFY: what does each order mean (order 1 is -2*lam*(t-tpeak)*exp(-lam*(t-tpeak)^2),
   *       order 0 uses exp(-2*lam*(t-tpeak)^2))? Is it a derivative order of the Ricker wavelet?
   * @param[in] tpeak Time of the pulse center.
   * @return Function value at time_n. No amplitude normalization is applied.
   */
  float evaluateRicker(float const& time_n, float const& f0, int order, float const& tpeak) {
    float pulse = 0.0;
    if ((time_n <= -0.9 * tpeak) || (time_n >= 2.9 * tpeak)) {
      return pulse;
    }

    constexpr float pi = M_PI;
    float const lam = (f0 * pi) * (f0 * pi);

    switch (order) {
      case 4: {
        pulse = 4.0 * lam * lam *
                (3.0 - 12.0 * lam * (time_n - tpeak) * (time_n - tpeak) +
                 4.0 * lam * lam * (time_n - tpeak) * (time_n - tpeak) * (time_n - tpeak) * (time_n - tpeak)) *
                exp(-lam * (time_n - tpeak) * (time_n - tpeak));
      } break;
      case 3: {
        pulse = 4.0 * lam * lam * (time_n - tpeak) * (3.0 - 2.0 * lam * (time_n - tpeak) * (time_n - tpeak)) *
                exp(-lam * (time_n - tpeak) * (time_n - tpeak));
      } break;
      case 2: {
        pulse = 2.0 * lam * (2.0 * lam * (time_n - tpeak) * (time_n - tpeak) - 1.0) *
                exp(-lam * (time_n - tpeak) * (time_n - tpeak));
      } break;
      case 1: {
        pulse = -2.0 * lam * (time_n - tpeak) * exp(-lam * (time_n - tpeak) * (time_n - tpeak));
      } break;
      case 0: {
        pulse = -(time_n - tpeak) * exp(-2 * lam * (time_n - tpeak) * (time_n - tpeak));
      } break;
      default:
        std::cout << "This option is not supported yet, rickerOrder must be 0, 1 or 2" << std::endl;
        break;
    }

    return pulse;
  }

  /**
   * @brief Samples evaluateRicker() at times i * timeSample.
   *
   * @param[in] nSamples Number of samples.
   * @param[in] timeSample Time step between samples.
   * @param[in] f0 Central frequency.
   * @param[in] order Function selector, see evaluateRicker().
   * @param[in] tpeak Time of the pulse center.
   * @return Vector of size nSamples, entry i is the value at time i * timeSample (first sample at time 0).
   */
  std::vector<float> computeSourceTerm(const int nSamples, const float timeSample, const float f0, const int order,
                                       const float tpeak) {
    std::vector<float> sourceTerm(nSamples);
    for (int i = 0; i < nSamples; i++) {
      float time_n = i * timeSample;
      sourceTerm[i] = evaluateRicker(time_n, f0, order, tpeak);
    }
    return sourceTerm;
  }
};
#endif  // FUNTIDES_ACQUISITION_INCLUDE_SOURCE_TIME_FUNCTION_H_
