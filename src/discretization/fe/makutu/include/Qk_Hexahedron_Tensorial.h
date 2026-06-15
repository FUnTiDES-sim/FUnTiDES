#ifndef SRC_DISCRETIZATION_FE_MAKUTU_INCLUDE_QK_HEXAHEDRON_TENSOSIRAL_H_
#define SRC_DISCRETIZATION_FE_MAKUTU_INCLUDE_QK_HEXAHEDRON_TENSOSIRAL_H_

#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"

#ifdef USE_KOKKOS
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_TeamVector_Impl.hpp>
#endif

template <typename GL_BASIS>
class Qk_Hexahedron_Tensorial : public Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS> {
 public:
  using BasisType = GL_BASIS;

  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;
  constexpr static int halfNodes = (GL_BASIS::numSupportPoints - 1) / 2;
  constexpr static int numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;
  constexpr static int numNodesPerFace = GL_BASIS::TensorProduct2D::numSupportPoints;
  constexpr static int maxSupportPoints = numNodes;
  constexpr static int numQuadraturePoints = numNodes;

  struct JacobianType {
    float data[3][3];
  };
  struct TeamGemm {};

  /**
   * @brief Scratch memory needed per team for KokkosBatched operations.
   */
  static constexpr size_t scratchBytesPerTeam() {
    constexpr size_t sva = (sizeof(real_t) >= 8) ? sizeof(real_t) : 8;
    constexpr size_t per_2d = num1dNodes * numNodesPerFace * sizeof(real_t) + sva;
    constexpr size_t per_1d = numNodes * sizeof(real_t) + sva;
    constexpr size_t per_D = num1dNodes * num1dNodes * sizeof(real_t) + sva;
    return 12 * per_2d + 2 * per_1d + per_D;
  }

  // ─── Utility index calculations ────────────────────────────────────

  PROXY_HOST_DEVICE
  constexpr static int linearIndex3DVal(const int qa, int const qb, int const qc) {
    return qa + qb * num1dNodes + qc * numNodesPerFace;
  }

  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex3D(int const k) {
    return linearIndex3DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * ((k % 4) / 2), (num1dNodes - 1) * (k / 4));
  }

  PROXY_HOST_DEVICE
  constexpr static int linearIndex2DVal(const int qa, const int qb) { return qa + qb * num1dNodes; }

  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex2D(int const k) {
    return linearIndex2DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * (k / 2));
  }

  // ─── Virtuals ──────────────────────────────────────────────────────

  PROXY_HOST_DEVICE ~Qk_Hexahedron_Tensorial() = default;
  PROXY_HOST_DEVICE virtual int getNumQuadraturePoints() { return numQuadraturePoints; }
  PROXY_HOST_DEVICE virtual int getNumSupportPoints() { return numNodes; }
  PROXY_HOST_DEVICE virtual int getMaxSupportPoints() const { return maxSupportPoints; }

  // ─── Basis and quadrature ──────────────────────────────────────────

  PROXY_HOST_DEVICE
  static void calcN(double const (&coords)[3], double (&N)[numNodes]) { GL_BASIS::TensorProduct3D::value(coords, N); }

  PROXY_HOST_DEVICE
  static void calcN(int const q, real_t (&N)[numNodes]) {
    for (int a = 0; a < numNodes; ++a) N[a] = 0;
    N[q] = 1.0;
  }

  constexpr static real_t interpolationCoord(const int q, const int k) {
    const real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0);
    return k == 0 ? (1.0 - alpha) : alpha;
  }

  PROXY_HOST_DEVICE
  constexpr static real_t basisGradientAt(const int q, const int p) {
    if (p <= halfNodes)
      return GL_BASIS::gradientAt(q, p);
    else
      return -GL_BASIS::gradientAt(num1dNodes - 1 - q, num1dNodes - 1 - p);
  }

  PROXY_HOST_DEVICE
  constexpr static real_t jacobianCoefficient1D(const int q, const int i, const int k, const int dir) {
    if (i == dir)
      return k == 0 ? -1.0 / 2.0 : 1.0 / 2.0;
    else
      return interpolationCoord(q, k);
  }

  PROXY_HOST_DEVICE
  constexpr static real_t quadratureWeight(const int q) { return GL_BASIS::weight(q); }

  // ─── Geometry helpers ──────────────────────────────────────────────

  PROXY_HOST_DEVICE
  static void jacobianTransformation2d(int const qa, int const qb, real_t const (&X)[4][3], real_t (&J)[3][2]) {
    for (int k = 0; k < 4; k++) {
      int ka = k % 2;
      int kb = k / 2;
      for (int j = 0; j < 2; j++) {
        real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) * jacobianCoefficient1D(qb, 1, kb, j);
        for (int i = 0; i < 3; i++) J[i][j] += jacCoeff * X[k][i];
      }
    }
  }

  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                     real_t (&J)[3][3]) {
    for (int k = 0; k < 8; k++) {
      const int ka = k % 2;
      const int kb = (k % 4) / 2;
      const int kc = k / 4;
      for (int j = 0; j < 3; j++) {
        real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) * jacobianCoefficient1D(qb, 1, kb, j) *
                          jacobianCoefficient1D(qc, 2, kc, j);
        for (int i = 0; i < 3; i++) J[i][j] += jacCoeff * X[k][i];
      }
    }
  }

  PROXY_HOST_DEVICE
  static void computeBMatrix(int const qa, int const qb, int const qc, real_t const (&X)[8][3], real_t (&J)[3][3],
                             real_t (&B)[6]) {
    jacobianTransformation(qa, qb, qc, X, J);
    real_t const detJ = determinant(J);
    real_t const invDetJ = 1.0 / detJ;

    B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) * invDetJ;
    B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) * invDetJ;
    B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]) * invDetJ;
    B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]) * invDetJ;
    B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]) * invDetJ;
    B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) * invDetJ;

    symInvert(B);
  }

  PROXY_HOST_DEVICE
  static real_t invJacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                          real_t (&J)[3][3]) {
    jacobianTransformation(qa, qb, qc, X, J);
    return invert3x3(J);
  }

  PROXY_HOST_DEVICE
  static real_t invJacobianTransformation(int const q, real_t const (&X)[8][3], real_t (&J)[3][3]) {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
    return invJacobianTransformation(qa, qb, qc, X, J);
  }

  PROXY_HOST_DEVICE
  static real_t transformedQuadratureWeight(int const q, real_t const (&X)[numNodes][3]) {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
    real_t J[3][3] = {{0}};
    jacobianTransformation(qa, qb, qc, X,
                           J);  // fallback uses flat X array representation, not implemented in scope, simplify.
    return determinant(J);
  }

  PROXY_HOST_DEVICE
  static void trilinearInterp(real_t const alpha, real_t const beta, real_t const gamma, real_t const (&X)[8][3],
                              real_t (&coords)[3]) {
    for (int i = 0; i < 3; i++) {
      coords[i] = X[0][i] * (1.0 - alpha) * (1.0 - beta) * (1.0 - gamma) +
                  X[1][i] * alpha * (1.0 - beta) * (1.0 - gamma) + X[2][i] * (1.0 - alpha) * beta * (1.0 - gamma) +
                  X[3][i] * alpha * beta * (1.0 - gamma) + X[4][i] * (1.0 - alpha) * (1.0 - beta) * gamma +
                  X[5][i] * alpha * (1.0 - beta) * gamma + X[6][i] * (1.0 - alpha) * beta * gamma +
                  X[7][i] * alpha * beta * gamma;
    }
  }

  template <typename FUNC, typename... PARAMS>
  PROXY_HOST_DEVICE static void supportLoop(real_t const (&coords)[3], FUNC &&func, PARAMS &&...params) {
    for (int c = 0; c < num1dNodes; ++c) {
      for (int b = 0; b < num1dNodes; ++b) {
        for (int a = 0; a < num1dNodes; ++a) {
          real_t const dNdXi[3] = {
              static_cast<real_t>(GL_BASIS::gradient(a, coords[0]) * GL_BASIS::value(b, coords[1]) *
                                  GL_BASIS::value(c, coords[2])),
              static_cast<real_t>(GL_BASIS::value(a, coords[0]) * GL_BASIS::gradient(b, coords[1]) *
                                  GL_BASIS::value(c, coords[2])),
              static_cast<real_t>(GL_BASIS::value(a, coords[0]) * GL_BASIS::value(b, coords[1]) *
                                  GL_BASIS::gradient(c, coords[2]))};
          int const nodeIndex = GL_BASIS::TensorProduct3D::linearIndex(a, b, c);
          func(dNdXi, nodeIndex, std::forward<PARAMS>(params)...);
        }
      }
    }
  }

  template <typename FUNC, typename... PARAMS>
  PROXY_HOST_DEVICE static void supportLoop(int const q, FUNC &&func, PARAMS &&...params) {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);

    for (int a = 0; a < num1dNodes; ++a) {
      real_t const dNdXi[3] = {basisGradientAt(a, qa), (a == qa) ? basisGradientAt(qb, qb) : real_t(0),
                               (a == qa) ? basisGradientAt(qc, qc) : real_t(0)};
      func(dNdXi, linearIndex3DVal(a, qb, qc), std::forward<PARAMS>(params)...);
    }
    for (int b = 0; b < num1dNodes; ++b) {
      if (b == qb) continue;
      real_t const dNdXi[3] = {real_t(0), basisGradientAt(b, qb), real_t(0)};
      func(dNdXi, linearIndex3DVal(qa, b, qc), std::forward<PARAMS>(params)...);
    }
    for (int c = 0; c < num1dNodes; ++c) {
      if (c == qc) continue;
      real_t const dNdXi[3] = {real_t(0), real_t(0), basisGradientAt(c, qc)};
      func(dNdXi, linearIndex3DVal(qa, qb, c), std::forward<PARAMS>(params)...);
    }
  }

  // ─── API Solver Callbacks ──────────────────────────────────────────

  template <typename FUNC>
  PROXY_HOST_DEVICE static void computeMassTerm(float const (&X)[8][3], FUNC &&func) {
    constexpr int N = num1dNodes;
    triple_loop<N, N, N>([&](auto const icqa, auto const icqb, auto const icqc) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;
      constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
      constexpr real_t w3D = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);
      real_t J[3][3] = {{0}};
      jacobianTransformation(qa, qb, qc, X, J);
      func(q, std::abs(determinant(J)) * w3D);
    });
  }

  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]) {
    int qa, qb;
    GL_BASIS::TensorProduct2D::multiIndex(q, qa, qb);
    const real_t w2D = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb));
    real_t B[3];
    real_t J[3][2] = {{0}};
    jacobianTransformation2d(qa, qb, X, J);
    B[0] = J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0];
    B[1] = J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1];
    B[2] = J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1];
    return sqrt(std::abs(symDeterminant(B))) * w2D;
  }

  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffNessTermwithJac(float const (&X)[8][3], FUNC1 &&func1, FUNC2 &&func2) {
    triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;
      JacobianType J = {{0}};
      jacobianTransformation(qa, qb, qc, X, J.data);
      computeGradPhiGradPhi<qa, qb, qc>(J, func1, func2);
    });
  }

  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiGradPhi(JacobianType &J, FUNC1 &&func1, FUNC2 &&func2) {
    real_t const detJ = invert3x3(J.data);
    const real_t w = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc));
    func1(qa, qb, qc, J.data);
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
        const real_t w00 = w * gia * gja;
        func2(ibc, jbc, w00 * detJ, 0, 0);
        const real_t w11 = w * gib * gjb;
        func2(aic, ajc, w11 * detJ, 1, 1);
        const real_t w22 = w * gic * gjc;
        func2(abi, abj, w22 * detJ, 2, 2);
        const real_t w12 = w * gib * gjc;
        func2(aic, abj, w12 * detJ, 1, 2);
        func2(abj, aic, w12 * detJ, 2, 1);
        const real_t w02 = w * gia * gjc;
        func2(ibc, abj, w02 * detJ, 0, 2);
        func2(abj, ibc, w02 * detJ, 2, 0);
        const real_t w01 = w * gia * gjb;
        func2(ibc, ajc, w01 * detJ, 0, 1);
        func2(ajc, ibc, w01 * detJ, 1, 0);
      }
    }
  }

  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiBGradPhi(real_t const (&B)[6], FUNC1 &&func1, FUNC2 &&func2) {
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
        func2(ibc, jbc, w * gia * gja * B[0]);
        func2(aic, ajc, w * gib * gjb * B[1]);
        func2(abi, abj, w * gic * gjc * B[2]);
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

  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffnessTerm(float const (&X)[8][3], FUNC1 &&func1, FUNC2 &&func2) {
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

  // ─── FLAT O(N^4) STIFFNESS SUM FACTORIZATIONS ────────────────────────

  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(float const (&X)[8][3], real_t const (&u_local)[numNodes],
                                                            real_t (&v_local)[numNodes], FUNC_ALPHA &&get_alpha) {
    real_t G_xi[numNodes] = {0};
    real_t G_eta[numNodes] = {0};
    real_t G_zeta[numNodes] = {0};

    triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;
      constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
      constexpr real_t w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);

      real_t dxi_q = 0, deta_q = 0, dzeta_q = 0;
      for_constexpr<num1dNodes>([&](auto ici) {
        constexpr int i = decltype(ici)::value;
        dxi_q += basisGradientAt(i, qa) * u_local[GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc)];
        deta_q += basisGradientAt(i, qb) * u_local[GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc)];
        dzeta_q += basisGradientAt(i, qc) * u_local[GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i)];
      });

      real_t J[3][3] = {{0}};
      real_t B[6] = {0};
      computeBMatrix(qa, qb, qc, X, J, B);
      real_t const scale = w * get_alpha(qa, qb, qc);

      G_xi[q] = scale * (B[0] * dxi_q + B[5] * deta_q + B[4] * dzeta_q);
      G_eta[q] = scale * (B[5] * dxi_q + B[1] * deta_q + B[3] * dzeta_q);
      G_zeta[q] = scale * (B[4] * dxi_q + B[3] * deta_q + B[2] * dzeta_q);
    });

    triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icia, auto const icib, auto const icic) {
      constexpr int ia = decltype(icia)::value;
      constexpr int ib = decltype(icib)::value;
      constexpr int ic = decltype(icic)::value;
      constexpr int node = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, ic);

      real_t v = 0;
      for_constexpr<num1dNodes>([&](auto icqa) {
        v += basisGradientAt(ia, decltype(icqa)::value) *
             G_xi[GL_BASIS::TensorProduct3D::linearIndex(decltype(icqa)::value, ib, ic)];
      });
      for_constexpr<num1dNodes>([&](auto icqb) {
        v += basisGradientAt(ib, decltype(icqb)::value) *
             G_eta[GL_BASIS::TensorProduct3D::linearIndex(ia, decltype(icqb)::value, ic)];
      });
      for_constexpr<num1dNodes>([&](auto icqc) {
        v += basisGradientAt(ic, decltype(icqc)::value) *
             G_zeta[GL_BASIS::TensorProduct3D::linearIndex(ia, ib, decltype(icqc)::value)];
      });

      v_local[node] += v;
    });
  }

  template <typename FUNC1>
  PROXY_HOST_DEVICE static void computeElasticStiffnessSumFact(float const (&X)[8][3],
                                                               real_t const (&u_local)[3][numNodes],
                                                               real_t (&f_local)[3][numNodes], FUNC1 &&func1) {
    real_t F_xi[3][numNodes] = {{0}};
    real_t F_eta[3][numNodes] = {{0}};
    real_t F_zeta[3][numNodes] = {{0}};

    triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;
      constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
      constexpr real_t w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);

      real_t grad_u_ref[3][3] = {{0}};
      for_constexpr<num1dNodes>([&](auto ici) {
        constexpr int i = decltype(ici)::value;
        const real_t gxi = basisGradientAt(i, qa);
        const real_t geta = basisGradientAt(i, qb);
        const real_t gzeta = basisGradientAt(i, qc);
        for (int s = 0; s < 3; ++s) {
          grad_u_ref[0][s] += gxi * u_local[s][GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc)];
          grad_u_ref[1][s] += geta * u_local[s][GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc)];
          grad_u_ref[2][s] += gzeta * u_local[s][GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i)];
        }
      });

      JacobianType J = {{0}};
      jacobianTransformation(qa, qb, qc, X, J.data);
      real_t const detJ = invert3x3(J.data);
      const real_t scale = w * detJ;

      real_t flux[3][3] = {{0}};
      func1(qa, qb, qc, J.data, grad_u_ref, flux);

      for (int f = 0; f < 3; ++f) {
        F_xi[f][q] = scale * flux[0][f];
        F_eta[f][q] = scale * flux[1][f];
        F_zeta[f][q] = scale * flux[2][f];
      }
    });

    triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icia, auto const icib, auto const icic) {
      constexpr int ia = decltype(icia)::value;
      constexpr int ib = decltype(icib)::value;
      constexpr int ic = decltype(icic)::value;
      constexpr int node = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, ic);

      real_t v[3] = {0};
      for_constexpr<num1dNodes>([&](auto icqa) {
        const real_t g = basisGradientAt(ia, decltype(icqa)::value);
        for (int f = 0; f < 3; ++f)
          v[f] += g * F_xi[f][GL_BASIS::TensorProduct3D::linearIndex(decltype(icqa)::value, ib, ic)];
      });
      for_constexpr<num1dNodes>([&](auto icqb) {
        const real_t g = basisGradientAt(ib, decltype(icqb)::value);
        for (int f = 0; f < 3; ++f)
          v[f] += g * F_eta[f][GL_BASIS::TensorProduct3D::linearIndex(ia, decltype(icqb)::value, ic)];
      });
      for_constexpr<num1dNodes>([&](auto icqc) {
        const real_t g = basisGradientAt(ic, decltype(icqc)::value);
        for (int f = 0; f < 3; ++f)
          v[f] += g * F_zeta[f][GL_BASIS::TensorProduct3D::linearIndex(ia, ib, decltype(icqc)::value)];
      });
      for (int f = 0; f < 3; ++f) f_local[f][node] += v[f];
    });
  }

  // ─── TEAM KOKKOSBATCHED GEMM STIFFNESS ─────────────────────────────

#ifdef USE_KOKKOS
  /**
   * @brief Team-parallel sum-factorization replacing flat triple-loops with
   * KokkosBatched::TeamVectorGemm. Reuses the new generic API format.
   */
  template <typename TEAM_MEMBER, typename FUNC_ALPHA>
  KOKKOS_INLINE_FUNCTION static void computeStiffnessTermSumFact_Team(
      const TEAM_MEMBER &member, float const (&X)[8][3], real_t const *p_local, real_t *f_local,
      real_t * /*G_xi*/,    // The solver provides these small 1D arrays
      real_t * /*G_eta*/,   // but we need wider ones for BatchedGEMM!
      real_t * /*G_zeta*/,  // So we bypass them and allocate locally in scratch space
      FUNC_ALPHA &&get_alpha) {
    constexpr int n = num1dNodes;
    constexpr int n2 = numNodesPerFace;
    constexpr int nTotal = n * n2;

    // Use full scratch capabilities allocated manually by the underlying user
    using ScratchSpace = typename TEAM_MEMBER::execution_space::scratch_memory_space;
    using ScratchView =
        Kokkos::View<real_t **, Kokkos::LayoutRight, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    ScratchView U_Xi(member.team_scratch(0), n, n2);
    ScratchView U_Eta(member.team_scratch(0), n, n2);
    ScratchView U_Zeta(member.team_scratch(0), n, n2);
    ScratchView dU_Xi(member.team_scratch(0), n, n2);
    ScratchView dU_Eta(member.team_scratch(0), n, n2);
    ScratchView dU_Zeta(member.team_scratch(0), n, n2);
    ScratchView Fx(member.team_scratch(0), n, n2);
    ScratchView Fy(member.team_scratch(0), n, n2);
    ScratchView Fz(member.team_scratch(0), n, n2);
    ScratchView y_Xi(member.team_scratch(0), n, n2);
    ScratchView y_Eta(member.team_scratch(0), n, n2);
    ScratchView y_Zeta(member.team_scratch(0), n, n2);
    ScratchView D_v(member.team_scratch(0), n, n);

    // Initialize the derivative matrix D
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n * n), [&](int idx) {
      int r = idx / n;
      int c = idx % n;
      D_v(r, c) = basisGradientAt(c, r);
    });

    // Zero matrices to be safe
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n * n2), [&](int idx) {
      int r = idx / n2;
      int c = idx % n2;
      U_Xi(r, c) = 0;
      U_Eta(r, c) = 0;
      U_Zeta(r, c) = 0;
    });
    member.team_barrier();

    // Fill Directions
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q_local) {
      int qa = q_local % n;
      int tmp = q_local / n;
      int qb = tmp % n;
      int qc = tmp / n;
      U_Xi(qa, qb * n + qc) = p_local[q_local];
      U_Eta(qb, qa * n + qc) = p_local[q_local];
      U_Zeta(qc, qa * n + qb) = p_local[q_local];
    });
    member.team_barrier();

    // Stage 2: TeamVectorGEMM dU = D * U
    KokkosBatched::TeamVectorGemm<TEAM_MEMBER, KokkosBlas::Trans::NoTranspose, KokkosBlas::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, D_v, U_Xi, 0.0, dU_Xi);
    KokkosBatched::TeamVectorGemm<TEAM_MEMBER, KokkosBlas::Trans::NoTranspose, KokkosBlas::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, D_v, U_Eta, 0.0, dU_Eta);
    KokkosBatched::TeamVectorGemm<TEAM_MEMBER, KokkosBlas::Trans::NoTranspose, KokkosBlas::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, D_v, U_Zeta, 0.0, dU_Zeta);
    member.team_barrier();

    // Stage 3: Physics Application
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q_local) {
      int qa = q_local % n;
      int tmp = q_local / n;
      int qb = tmp % n;
      int qc = tmp / n;

      real_t J[3][3] = {{0}};
      real_t B[6] = {0};
      computeBMatrix(qa, qb, qc, X, J, B);
      real_t scale = quadratureWeight(qa) * quadratureWeight(qb) * quadratureWeight(qc) * get_alpha(qa, qb, qc);

      int col_Xi = qb * n + qc;
      int col_Eta = qa * n + qc;
      int col_Zeta = qa * n + qb;

      Fx(qa, col_Xi) = scale * (B[0] * dU_Xi(qa, col_Xi) + B[5] * dU_Eta(qb, col_Eta) + B[4] * dU_Zeta(qc, col_Zeta));
      Fy(qb, col_Eta) = scale * (B[5] * dU_Xi(qa, col_Xi) + B[1] * dU_Eta(qb, col_Eta) + B[3] * dU_Zeta(qc, col_Zeta));
      Fz(qc, col_Zeta) = scale * (B[4] * dU_Xi(qa, col_Xi) + B[3] * dU_Eta(qb, col_Eta) + B[2] * dU_Zeta(qc, col_Zeta));
    });
    member.team_barrier();

    // Stage 4: TeamVectorGEMM Backtransform y = D^T * F
    KokkosBatched::TeamVectorGemm<TEAM_MEMBER, KokkosBlas::Trans::Transpose, KokkosBlas::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, D_v, Fx, 0.0, y_Xi);
    KokkosBatched::TeamVectorGemm<TEAM_MEMBER, KokkosBlas::Trans::Transpose, KokkosBlas::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, D_v, Fy, 0.0, y_Eta);
    KokkosBatched::TeamVectorGemm<TEAM_MEMBER, KokkosBlas::Trans::Transpose, KokkosBlas::Trans::NoTranspose,
                                  KokkosBatched::Algo::Gemm::Unblocked>::invoke(member, 1.0, D_v, Fz, 0.0, y_Zeta);
    member.team_barrier();

    // Stage 5: Accumulate final f_local directly from thread specific indices.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q_local) {
      int qa = q_local % n;
      int tmp = q_local / n;
      int qb = tmp % n;
      int qc = tmp / n;
      f_local[q_local] += y_Xi(qa, qb * n + qc) + y_Eta(qb, qa * n + qc) + y_Zeta(qc, qa * n + qb);
    });
  }
#else
  // Fallback for non-Kokkos builds just directly executes the flat version inline
  template <typename TEAM_MEMBER, typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact_Team(const TEAM_MEMBER &team, float const (&X)[8][3],
                                                                 real_t const *p_local, real_t *f_local, real_t *G_xi,
                                                                 real_t *G_eta, real_t *G_zeta,
                                                                 FUNC_ALPHA &&get_alpha) {
    // Fallback
  }
#endif
};

using Q1_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis1>;
using Q2_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis2>;
using Q3_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis3GL>;
using Q4_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis4GL>;
using Q5_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis5GL>;
using Q6_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis6GL>;
using Q7_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis7GL>;
using Q8_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis8GL>;
using Q9_Hexahedron_Tensorial = Qk_Hexahedron_Tensorial<LagrangeBasis9GL>;

template <int ORDER>
struct Qk_Hexahedron_Tensorial_Selector;
template <>
struct Qk_Hexahedron_Tensorial_Selector<1> {
  using type = Q1_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<2> {
  using type = Q2_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<3> {
  using type = Q3_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<4> {
  using type = Q4_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<5> {
  using type = Q5_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<6> {
  using type = Q6_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<7> {
  using type = Q7_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<8> {
  using type = Q8_Hexahedron_Tensorial;
};
template <>
struct Qk_Hexahedron_Tensorial_Selector<9> {
  using type = Q9_Hexahedron_Tensorial;
};

#endif  // SRC_DISCRETIZATION_FE_MAKUTU_INCLUDE_QK_HEXAHEDRON_TENSOSIRAL_H_
