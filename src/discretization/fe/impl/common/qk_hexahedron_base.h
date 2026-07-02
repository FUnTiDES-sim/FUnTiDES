#ifndef FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_BASE_H_
#define FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_BASE_H_

#include <data_type.h>

#include "fe_discretization.h"
#include "mathUtilites.h"

/**
 * @brief Implementation-independent core of the Qk hexahedral Gauss-Lobatto
 *        spectral-element discretization.
 *
 * Holds every quantity and kernel shared by all back-ends: compile-time sizes,
 * index maps, the isoparametric Jacobian / B-matrix, and the mass, damping and
 * acoustic-stiffness (flat sum-factorization) kernels. The concrete back-ends
 * (makutu Lagrange-GaussLobatto and tensorial GEMM) derive from this class and
 * add only their own machinery.
 *
 * Every member is static and host/device-callable: no vtable, no instance
 * state. This is mandatory for GPU execution -- all dispatch is compile-time.
 */
template <typename GL_BASIS>
class QkHexahedronBase : public discretization::fe::api::FeDiscretizationTag {
 public:
  // Expose the basis type for tests and external use
  using BasisType = GL_BASIS;

  /// The number of nodes/support points per element per dimension.
  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;

  /// Half the number of support points, rounded down. Precomputed for
  /// efficiency
  constexpr static int halfNodes = (GL_BASIS::numSupportPoints - 1) / 2;

  /// The number of nodes/support points per element.
  constexpr static int numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;

  /// The number of nodes/support points per face
  constexpr static int numNodesPerFace = GL_BASIS::TensorProduct2D::numSupportPoints;

  /// The maximum number of support points per element.
  constexpr static int maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static int numQuadraturePoints = numNodes;

  struct JacobianType {
    float data[3][3];
  };

  /**
   * @brief The linear index associated to the given one-dimensional indices in
   * the three directions
   * @param qa The index in the first direction
   * @param qb The index in the second direction
   * @param qc The index in the third direction
   * @return The linear index in 3D
   */
  PROXY_HOST_DEVICE
  constexpr static int linearIndex3DVal(const int qa, int const qb, int const qc) {
    return qa + qb * num1dNodes + qc * numNodesPerFace;
  }

  /**
   * @brief Converts from the index of the point in the mesh and the linear 3D
   * index of the corresponding dof.
   * @param k The index of the mesh vertex, from 0 to 7
   * @return The linear index in 3D
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex3D(int const k) {
    return linearIndex3DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * ((k % 4) / 2), (num1dNodes - 1) * (k / 4));
  }

  /**
   * @brief The linear index associated to the given one-dimensional indices in
   * the two directions
   * @param qa The index in the first direction
   * @param qb The index in the second direction
   * @return The linear index in 2D
   */
  PROXY_HOST_DEVICE
  constexpr static int linearIndex2DVal(const int qa, const int qb) { return qa + qb * num1dNodes; }

  /**
   * @brief Converts from the index of the point in the mesh and the linear 2D
   * index of the corresponding dof.
   * @param k The index of the mesh vertex, from 0 to 3
   * @return The linear index in 2D
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex2D(int const k) {
    return linearIndex2DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * (k / 2));
  }

  /**
   * @brief Compute the interpolation coefficients of the q-th quadrature point
   * in a given direction
   * @param q the index of the quadrature point in 1D
   * @param k the index of the interval endpoint (0 or 1)
   * @return The interpolation coefficient
   */
  constexpr static real_t interpolationCoord(const int q, const int k) {
    const real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0);
    return k == 0 ? (1.0 - alpha) : alpha;
  }

  /**
   * @brief Compute the 1st derivative of the q-th 1D basis function at
   * quadrature point p
   * @param q the index of the 1D basis funcion
   * @param p the index of the 1D quadrature point
   * @return The derivative value
   */
  PROXY_HOST_DEVICE
  constexpr static real_t basisGradientAt(const int q, const int p) {
    if (p <= halfNodes) {
      return GL_BASIS::gradientAt(q, p);
    } else {
      return -GL_BASIS::gradientAt(GL_BASIS::numSupportPoints - 1 - q, GL_BASIS::numSupportPoints - 1 - p);
    }
  }

  /**
   * @brief Compute the 1D factor of the coefficient of the jacobian on the q-th
   * quadrature point, with respect to the k-th interval endpoint (0 or 1). The
   * computation depends on the position in the basis tensor product of this
   * term (i, equal to 0, 1 or 2) and on the direction in which the gradient is
   * being computed (dir, from 0 to 2)
   * @param q The index of the quadrature point in 1D
   * @param i The index of the position in the tensor product
   * @param k The index of the interval endpoint (0 or 1)
   * @param dir The direction in which the derivatives are being computed
   * @return The value of the jacobian factor
   */
  PROXY_HOST_DEVICE
  constexpr static real_t jacobianCoefficient1D(const int q, const int i, const int k, const int dir) {
    if (i == dir) {
      return k == 0 ? -1.0 / 2.0 : 1.0 / 2.0;
    } else {
      return interpolationCoord(q, k);
    }
  }

  /// 1D Gauss-Lobatto quadrature weight at node q.
  PROXY_HOST_DEVICE
  constexpr static real_t quadratureWeight(const int q) { return GL_BASIS::weight(q); }

  // Non-virtual, compile-time replacements for the former virtual hooks.
  // (The old virtual versions added a vtable that was never dispatched and is
  //  unsafe on device.)
  PROXY_HOST_DEVICE static constexpr int getNumQuadraturePoints() { return numQuadraturePoints; }
  PROXY_HOST_DEVICE static constexpr int getNumSupportPoints() { return numNodes; }
  PROXY_HOST_DEVICE static constexpr int getMaxSupportPoints() { return maxSupportPoints; }

  PROXY_HOST_DEVICE
  static void jacobianTransformation2d(int const qa, int const qb, real_t const (&X)[4][3], real_t (&J)[3][2]);

  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                     real_t (&J)[3][3]);

  template <typename FUNC>
  PROXY_HOST_DEVICE static void computeMassTerm(float const (&X)[8][3], FUNC &&func);

  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]);

  PROXY_HOST_DEVICE
  static void computeBMatrix(int const qa, int const qb, int const qc, real_t const (&X)[8][3], real_t (&J)[3][3],
                             real_t (&B)[6]);

  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffnessTerm(float const (&X)[8][3], FUNC1 &&func1, FUNC2 &&func2);

  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(float const (&X)[8][3], real_t const (&p_local)[numNodes],
                                                            real_t (&f_local)[numNodes], FUNC_ALPHA &&get_alpha);

  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiBGradPhi(real_t const (&B)[6], FUNC1 &&func1, FUNC2 &&func2);

 protected:
  // Static-interface base: never deleted through a base pointer, never holds a
  // vtable. Protected + non-virtual is the GPU-safe choice.
  PROXY_HOST_DEVICE ~QkHexahedronBase() = default;
};

template <typename GL_BASIS>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::jacobianTransformation(int const qa, int const qb, int const qc,
                                                                          real_t const (&X)[8][3], real_t (&J)[3][3]) {
  for (int k = 0; k < 8; k++) {
    const int ka = k % 2;
    const int kb = (k % 4) / 2;
    const int kc = k / 4;
    for (int j = 0; j < 3; j++) {
      real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) * jacobianCoefficient1D(qb, 1, kb, j) *
                        jacobianCoefficient1D(qc, 2, kc, j);
      for (int i = 0; i < 3; i++) {
        J[i][j] += jacCoeff * X[k][i];
      }
    }
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::jacobianTransformation2d(int const qa, int const qb,
                                                                            real_t const (&X)[4][3],
                                                                            real_t (&J)[3][2]) {
  for (int k = 0; k < 4; k++) {
    int ka = k % 2;
    int kb = k / 2;
    for (int j = 0; j < 2; j++) {
      real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) * jacobianCoefficient1D(qb, 1, kb, j);
      for (int i = 0; i < 3; i++) {
        J[i][j] += jacCoeff * X[k][i];
      }
    }
  }
}

template <typename GL_BASIS>
template <typename FUNC>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::computeMassTerm(float const (&X)[8][3], FUNC &&func) {
  constexpr int N = num1dNodes;
  triple_loop<N, N, N>([&](auto const icqa, auto const icqb, auto const icqc) {
    constexpr int qa = decltype(icqa)::value;
    constexpr int qb = decltype(icqb)::value;
    constexpr int qc = decltype(icqc)::value;
    constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
    constexpr real_t w3D = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);
    real_t J[3][3] = {{0}};
    jacobianTransformation(qa, qb, qc, X, J);
    real_t val = std::abs(determinant(J)) * w3D;
    func(q, val);
  });
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t QkHexahedronBase<GL_BASIS>::computeDampingTerm(int const q, real_t const (&X)[4][3]) {
  int qa, qb;
  GL_BASIS::TensorProduct2D::multiIndex(q, qa, qb);
  const real_t w2D = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb));
  real_t B[3];
  real_t J[3][2] = {{0}};
  jacobianTransformation2d(qa, qb, X, J);
  // compute J^T.J, using Voigt notation for B
  B[0] = J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0];
  B[1] = J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1];
  B[2] = J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1];
  return sqrt(std::abs(symDeterminant(B))) * w2D;
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::computeBMatrix(int const qa, int const qb, int const qc,
                                                                  real_t const (&X)[8][3], real_t (&J)[3][3],
                                                                  real_t (&B)[6]) {
  jacobianTransformation(qa, qb, qc, X, J);
  real_t const detJ = determinant(J);
  real_t const invDetJ = 1.0 / detJ;

  // compute J^T.J/det(J), using Voigt notation for B
  B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) * invDetJ;
  B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) * invDetJ;
  B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]) * invDetJ;
  B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]) * invDetJ;
  B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]) * invDetJ;
  B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) * invDetJ;

  // compute detJ*J^{-1}J^{-T}
  symInvert(B);
}

template <typename GL_BASIS>
template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::computeGradPhiBGradPhi(real_t const (&B)[6], FUNC1 &&func1,
                                                                          FUNC2 &&func2) {
  const real_t w = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc));
  func1(qa, qb, qc);
  for (int i = 0; i < num1dNodes; i++) {
    const int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
    const int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
    const int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
    const real_t gia = basisGradientAt(i, qa);
    const real_t gib = basisGradientAt(i, qb);
    const real_t gic = basisGradientAt(i, qc);
    for (int j = 0; j < num1dNodes; j++) {
      const int jbc = GL_BASIS::TensorProduct3D::linearIndex(j, qb, qc);
      const int ajc = GL_BASIS::TensorProduct3D::linearIndex(qa, j, qc);
      const int abj = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, j);
      const real_t gja = basisGradientAt(j, qa);
      const real_t gjb = basisGradientAt(j, qb);
      const real_t gjc = basisGradientAt(j, qc);
      // diagonal terms
      const real_t w0 = w * gia * gja;
      func2(ibc, jbc, w0 * B[0]);
      const real_t w1 = w * gib * gjb;
      func2(aic, ajc, w1 * B[1]);
      const real_t w2 = w * gic * gjc;
      func2(abi, abj, w2 * B[2]);
      // off-diagonal terms
      const real_t w3 = w * gib * gjc;
      func2(aic, abj, w3 * B[3]);
      func2(abj, aic, w3 * B[3]);
      const real_t w4 = w * gia * gjc;
      func2(ibc, abj, w4 * B[4]);
      func2(abj, ibc, w4 * B[4]);
      const real_t w5 = w * gia * gjb;
      func2(ibc, ajc, w5 * B[5]);
      func2(ajc, ibc, w5 * B[5]);
    }
  }
}

template <typename GL_BASIS>
template <typename FUNC1, typename FUNC2>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::computeStiffnessTerm(float const (&X)[8][3], FUNC1 &&func1,
                                                                        FUNC2 &&func2) {
  triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
    constexpr int qa = decltype(icqa)::value;
    constexpr int qb = decltype(icqb)::value;
    constexpr int qc = decltype(icqc)::value;
    real_t B[6] = {0};
    real_t J[3][3] = {{0}};
    computeBMatrix(qa, qb, qc, X, J, B);
    computeGradPhiBGradPhi<qa, qb, qc>(B, func1, func2);
  });
}

template <typename GL_BASIS>
template <typename FUNC_ALPHA>
PROXY_HOST_DEVICE void QkHexahedronBase<GL_BASIS>::computeStiffnessTermSumFact(float const (&X)[8][3],
                                                                               real_t const (&u_local)[numNodes],
                                                                               real_t (&v_local)[numNodes],
                                                                               FUNC_ALPHA &&get_alpha) {
  // Weighted fluxes G^{ξ,η,ζ}[q] = w_q * alpha_q * M(B_q) · ∇_ξ u_q
  real_t G_xi[numNodes] = {0};
  real_t G_eta[numNodes] = {0};
  real_t G_zeta[numNodes] = {0};

  // Pass 1+2 fused: gradient of u, then application of the metric, alpha and weight
  triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
    constexpr int qa = decltype(icqa)::value;
    constexpr int qb = decltype(icqb)::value;
    constexpr int qc = decltype(icqc)::value;
    constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);

    // Quadrature weight management is internal to the math library
    constexpr real_t w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);

    real_t dxi_q = 0, deta_q = 0, dzeta_q = 0;
    for_constexpr<num1dNodes>([&](auto ici) {
      constexpr int i = decltype(ici)::value;
      constexpr int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
      constexpr int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
      constexpr int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);

      dxi_q += basisGradientAt(i, qa) * u_local[ibc];
      deta_q += basisGradientAt(i, qb) * u_local[aic];
      dzeta_q += basisGradientAt(i, qc) * u_local[abi];
    });

    real_t J[3][3] = {{0}};
    real_t B[6] = {0};
    computeBMatrix(qa, qb, qc, X, J, B);

    // 'scale' combines physics (alpha) and quadrature weight (w)
    real_t const scale = w * get_alpha(qa, qb, qc);

    G_xi[q] = scale * (B[0] * dxi_q + B[5] * deta_q + B[4] * dzeta_q);
    G_eta[q] = scale * (B[5] * dxi_q + B[1] * deta_q + B[3] * dzeta_q);
    G_zeta[q] = scale * (B[4] * dxi_q + B[3] * deta_q + B[2] * dzeta_q);
  });

  // Pass 3: divergence — v_{ia,ib,ic} += D^T·G^ξ + D^T·G^η + D^T·G^ζ
  triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icia, auto const icib, auto const icic) {
    constexpr int ia = decltype(icia)::value;
    constexpr int ib = decltype(icib)::value;
    constexpr int ic = decltype(icic)::value;
    constexpr int node = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, ic);

    real_t v = 0;
    for_constexpr<num1dNodes>([&](auto icqa) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int q_xi = GL_BASIS::TensorProduct3D::linearIndex(qa, ib, ic);
      v += basisGradientAt(ia, qa) * G_xi[q_xi];
    });
    for_constexpr<num1dNodes>([&](auto icqb) {
      constexpr int qb = decltype(icqb)::value;
      constexpr int q_eta = GL_BASIS::TensorProduct3D::linearIndex(ia, qb, ic);
      v += basisGradientAt(ib, qb) * G_eta[q_eta];
    });
    for_constexpr<num1dNodes>([&](auto icqc) {
      constexpr int qc = decltype(icqc)::value;
      constexpr int q_zeta = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, qc);
      v += basisGradientAt(ic, qc) * G_zeta[q_zeta];
    });

    v_local[node] += v;
  });
}

#endif  // FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_BASE_H_
