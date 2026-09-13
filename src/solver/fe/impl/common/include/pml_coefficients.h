#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_PML_COEFFICIENTS_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_PML_COEFFICIENTS_H_

#include <cmath>

namespace solver {
namespace fe {

/**
 * @brief Per-node Convolutional PML (C-PML) coefficients for the second-order
 *        acoustic wave equation.
 *
 * The acoustic pressure equation
 *
 *   1/kappa * d2p/dt2 = div( 1/rho * grad p )
 *
 * is solved in the PML layer with stretched spatial derivatives
 *
 *   d/dx_i  ->  (1/s_i) d/dx_i,      s_i = kappa_i + d_i / (alpha_i + i*omega)
 *
 * Following Komatitsch & Martin (2007) the stretched gradient of the trial
 * field is
 *
 *   d~p/dx_i = (1/kappa_i) dp/dx_i - (1/kappa_i) psi_i
 *
 * where the memory variable psi_i satisfies
 *
 *   d psi_i/dt + (alpha_i + d_i/kappa_i) psi_i = (d_i/kappa_i) dp/dx_i
 *
 * and is advanced with the exact convolution (Wang, Lee & Teixeira 2006,
 * eq. 21) using the current and previous gradient:
 *
 *   psi_i^{n+1} = coef0_i psi_i^n + coef1_i (dp/dx_i)^n + coef2_i (dp/dx_i)^{n-1}
 *
 * with
 *
 *   coef0_i = exp( -(alpha_i + d_i/kappa_i) dt )
 *   coef1_i = coef2_i = (d_i/kappa_i)/(alpha_i + d_i/kappa_i) * (1 - coef0_i)
 *
 * The profiles follow Komatitsch & Martin (2007): with delta the distance
 * from the PML inner boundary and L the layer thickness,
 *
 *   d_i(delta)     = d_max * (delta/L)^N
 *   kappa_i(delta) = 1 + (kappa_max - 1) * (delta/L)^N
 *   alpha_i(delta) = alpha_max * (1 - delta/L)^N
 *
 *   d_max = -(N+1) * vp / (2 L) * ln(R)
 *
 * where R is the target reflection coefficient and vp the local P velocity.
 *
 * With a zero profile (d = 0, kappa = 1) every coefficient reduces to the
 * identity: coef0 = 1, coef1 = coef2 = 0, so psi stays zero and the stretched
 * gradient equals the unstretched one.
 */
struct PmlCoefficients {
  // Stretching profile per direction (0=x, 1=y, 2=z).
  float d[3] = {0.0f, 0.0f, 0.0f};
  float kappa[3] = {1.0f, 1.0f, 1.0f};
  float alpha[3] = {0.0f, 0.0f, 0.0f};

  // Convolution coefficients per direction (Wang/Lee/Teixeira 2006).
  float coef0[3] = {1.0f, 1.0f, 1.0f};
  float coef1[3] = {0.0f, 0.0f, 0.0f};
  float coef2[3] = {0.0f, 0.0f, 0.0f};

  // True if this node lies inside the PML layer (any direction active).
  bool isPml = false;
};

/**
 * @brief Fill the C-PML coefficients for one node from its coordinates.
 *
 * The domain is assumed to span [0, domainSize[i]] in each direction. A node
 * is inside the layer of direction i when its coordinate is within pmlSize[i]
 * of either boundary; delta is then the distance from the inner edge of the
 * layer (0 at the inner edge, pmlSize[i] at the outer boundary).
 *
 * @param x, y, z    Node coordinates.
 * @param domainSize Domain extent in each direction (lx, ly, lz).
 * @param pmlSize    PML thickness in each direction (0 = no PML there).
 * @param dt         Time step.
 * @param vp         P-wave velocity at the node (scales d_max).
 * @param profile    Profile exponent N (default 2, quadratic).
 * @param reflection Target reflection coefficient R (default 1e-3).
 * @param alphaMax   Maximum alpha (default 0 = no frequency shift).
 * @param kappaMax   Maximum kappa (default 1 = no coordinate stretching).
 * @param out        Output coefficients.
 */
inline void fillPmlCoefficients(float x, float y, float z, const float domainSize[3], const float pmlSize[3],
                                float dt, float vp, float profile, float reflection, float alphaMax, float kappaMax,
                                PmlCoefficients& out) {
  const float coord[3] = {x, y, z};
  const float kProfile = (profile > 0.0f) ? profile : 2.0f;

  bool anyPml = false;
  for (int i = 0; i < 3; ++i) {
    const float L = pmlSize[i];
    const float D = domainSize[i];
    const float c = coord[i];

    // Distance from the inner PML boundary (0 = inner edge, L = outer edge).
    float delta = -1.0f;
    if (L > 0.0f && c < L) {
      delta = L - c;  // near the low-coordinate boundary
    } else if (L > 0.0f && c > D - L) {
      delta = c - (D - L);  // near the high-coordinate boundary
    }

    if (delta < 0.0f) {
      // Interior of the domain in this direction: identity coefficients.
      out.d[i] = 0.0f;
      out.kappa[i] = 1.0f;
      out.alpha[i] = 0.0f;
      out.coef0[i] = 1.0f;
      out.coef1[i] = 0.0f;
      out.coef2[i] = 0.0f;
      continue;
    }

    const float r = delta / L;  // 0 at inner edge, 1 at outer boundary
    const float rN = std::pow(r, kProfile);
    const float dMax = (vp > 0.0f && reflection > 0.0f)
                           ? -(kProfile + 1.0f) * vp / (2.0f * L) * std::log(reflection)
                           : 0.0f;

    out.d[i] = dMax * rN;
    out.kappa[i] = 1.0f + (kappaMax - 1.0f) * rN;
    out.alpha[i] = alphaMax * std::pow(1.0f - r, kProfile);

    // Convolution coefficients (Wang/Lee/Teixeira 2006, eq. 21).
    const float a = out.alpha[i] + out.d[i] / out.kappa[i];
    const float c0 = std::exp(-a * dt);
    out.coef0[i] = c0;
    if (a > 0.0f) {
      const float c1 = (out.d[i] / out.kappa[i]) / a * (1.0f - c0);
      out.coef1[i] = c1;
      out.coef2[i] = c1;
    } else {
      out.coef1[i] = 0.0f;
      out.coef2[i] = 0.0f;
    }

    anyPml = true;
  }

  out.isPml = anyPml;
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_PML_COEFFICIENTS_H_
