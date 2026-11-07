#pragma once

#include <algorithm>
#include <cmath>

/**
 * @file sponge_profile.h
 * @brief Damping profile of the absorbing sponge layers (shared by the SEM solvers).
 *
 * In a sponge of width L, both time levels of the field are multiplied at each step by
 * exp(-sigma(d) dt), with the damping rate
 *
 *     sigma(d) = sigma_max ((L - d) / L)^2,   0 <= d < L,
 *
 * d being the distance to the domain edge (sigma = 0 at the inner side of the sponge, with a
 * zero slope there, which limits the reflection on the sponge itself). sigma_max is set so that
 * a wave crossing the sponge at normal incidence at the velocity vmax and coming back is
 * attenuated by the factor R:
 *
 *     exp(-2 / vmax * integral_0^L sigma) = R   =>   sigma_max = 3 vmax ln(1 / R) / (2 L).
 *
 * Slower waves are attenuated more. Being a rate, sigma does not depend on the time step.
 * On a corner, the rates of the two (three) sponges add up.
 */
namespace sponge {

/// Attenuation used when the requested one is not in (0, 1).
inline constexpr double kDefaultReflection = 1e-2;

/**
 * @brief Damping rate (1/s) at the distance d (m) from the domain edge, in a sponge of width L (m).
 * @param reflection Amplitude ratio R after a normal round trip through the sponge at vmax.
 * @return 0 outside the sponge (d >= L) or when L <= 0.
 */
inline double dampingRate(double d, double width, double vmax, double reflection) {
  if (width <= 0.0 || d >= width) return 0.0;
  double const r = (reflection > 0.0 && reflection < 1.0) ? reflection : kDefaultReflection;
  double const s = (width - std::max(d, 0.0)) / width;
  return 1.5 * vmax * std::log(1.0 / r) / width * s * s;
}

}  // namespace sponge
