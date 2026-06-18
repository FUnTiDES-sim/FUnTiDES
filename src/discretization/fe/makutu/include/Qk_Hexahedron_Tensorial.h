#ifndef FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_TENSORIAL_GEMM_H_
#define FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_TENSORIAL_GEMM_H_

// data_type.h first: defines PROXY_HOST_DEVICE / real_t before the basis files.
#include <data_type.h>

#include "LagrangeBasis1.h"
#include "LagrangeBasis2.h"
#include "LagrangeBasis3GL.h"
#include "LagrangeBasis4GL.h"
#include "LagrangeBasis5GL.h"
#include "LagrangeBasis6GL.h"
#include "LagrangeBasis7GL.h"
#include "LagrangeBasis8GL.h"
#include "LagrangeBasis9GL.h"
#include "mathUtilites.h"  // determinant, invert3x3, symDeterminant, triple_loop, for_constexpr, ...

/**
 * @class Qk_Hexahedron_Tensorial_GEMM
 * @tparam GL_BASIS The Lagrange basis type (LagrangeBasis1, LagrangeBasis2, ...)
 */
template <typename GL_BASIS>
class Qk_Hexahedron_Tensorial_GEMM final {
 public:
  // Expose the basis type for tests and external use (makutu compatibility).
  using BasisType = GL_BASIS;

  /// Number of nodes/support points per element per dimension.
  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;

  /// Half the number of support points, rounded down (precomputed).
  constexpr static int halfNodes = (GL_BASIS::numSupportPoints - 1) / 2;

  /// Total number of nodes/support points per element.
  constexpr static int numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;

  /// Number of nodes/support points per face.
  constexpr static int numNodesPerFace = num1dNodes * num1dNodes;

  /// Maximum number of support points per element.
  constexpr static int maxSupportPoints = numNodes;

  /// Number of quadrature points per element (== nodes for GL).
  constexpr static int numQuadraturePoints = numNodes;

  struct JacobianType {
    float data[3][3];
  };

  /// Marker trait: an INTEGRAL_TYPE exposing a team GEMM stiffness path.
  struct TeamGemm {};

  //==========================================================================
  // Index helpers (ported verbatim from makutu — required by the solver and by
  // the assembly loops; identical node ordering: x fastest).
  //==========================================================================

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

  //==========================================================================
  // Basis / quadrature primitives
  //==========================================================================

  /// d(phi_p)/d(xi) evaluated at xi_q (direct polynomial eval, GL symmetry).
  PROXY_HOST_DEVICE
  constexpr static real_t basisGradientAt(const int q, const int p) {
    if (p <= halfNodes) {
      return GL_BASIS::gradientAt(q, p);
    } else {
      return -GL_BASIS::gradientAt(num1dNodes - 1 - q, num1dNodes - 1 - p);
    }
  }

  /// 1D Gauss-Lobatto quadrature weight at node q.
  PROXY_HOST_DEVICE
  constexpr static real_t quadratureWeight(const int q) { return GL_BASIS::weight(q); }

  PROXY_HOST_DEVICE
  constexpr static real_t interpolationCoord(const int q, const int k) {
    const real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0);
    return k == 0 ? (1.0 - alpha) : alpha;
  }

  PROXY_HOST_DEVICE
  constexpr static real_t jacobianCoefficient1D(const int q, const int i, const int k, const int dir) {
    if (i == dir)
      return k == 0 ? -0.5 : 0.5;
    else
      return interpolationCoord(q, k);
  }

  //==========================================================================
  // Jacobians
  //==========================================================================

  /// 3D isoparametric (trilinear, 8-vertex) Jacobian at GL point (qa,qb,qc).
  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                     real_t (&J)[3][3]) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) J[i][j] = 0.0;

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

  /// 2D face Jacobian (bilinear, 4-corner). Ported verbatim from makutu.
  PROXY_HOST_DEVICE
  static void jacobianTransformation2d(int const qa, int const qb, real_t const (&X)[4][3], real_t (&J)[3][2]) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 2; ++j) J[i][j] = 0.0;
    for (int k = 0; k < 4; k++) {
      int ka = k % 2;
      int kb = k / 2;
      for (int j = 0; j < 2; j++) {
        real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) * jacobianCoefficient1D(qb, 1, kb, j);
        for (int i = 0; i < 3; i++) J[i][j] += jacCoeff * X[k][i];
      }
    }
  }

  /**
   * @brief B = det(J) * (J^T J)^{-1} in Voigt notation [xx,yy,zz,yz,xz,xy].
   *
   * Same quantity and same Voigt ordering as makutu::computeBMatrix, so the
   * flux contraction below matches the makutu kernel term by term.
   */
  PROXY_HOST_DEVICE
  static void computeBMatrix(int const qa, int const qb, int const qc, real_t const (&X)[8][3], real_t (&J)[3][3],
                             real_t (&B)[6]) {
    jacobianTransformation(qa, qb, qc, X, J);
    real_t const detJ = determinant(J);
    real_t const invDetJ = 1.0 / detJ;

    // J^T.J / det(J), Voigt notation [xx,yy,zz,yz,xz,xy]
    B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) * invDetJ;
    B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) * invDetJ;
    B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]) * invDetJ;
    B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]) * invDetJ;
    B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]) * invDetJ;
    B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) * invDetJ;

    // B <- det(J) * (J^T J)^{-1}  (same quantity & Voigt order as makutu)
    symInvert(B);
  }

  //==========================================================================
  // MASS TERM  (ported verbatim from makutu — numerically identical)
  //==========================================================================

  /**
   * @brief Diagonal mass contribution per quadrature point.
   * @param X    8 corner coordinates of the element.
   * @param func Callback func(int q, real_t val).
   */
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
      real_t val = std::abs(determinant(J)) * w3D;
      func(q, val);
    });
  }

  //==========================================================================
  // DAMPING TERM  (ported verbatim from makutu — numerically identical)
  //==========================================================================

  /**
   * @brief Diagonal damping (face-integrated) contribution for d.o.f. q.
   * @param q Face quadrature point index (2D).
   * @param X 4 corner coordinates of the face.
   * @return  The surface quadrature factor sqrt(|det(J^T J)|) * w2D.
   */
  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]) {
    int qa, qb;
    GL_BASIS::TensorProduct2D::multiIndex(q, qa, qb);
    const real_t w2D = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb));
    real_t B[3];
    real_t J[3][2] = {{0}};
    jacobianTransformation2d(qa, qb, X, J);
    // J^T.J in Voigt notation (2x2 symmetric)
    B[0] = J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0];
    B[1] = J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1];
    B[2] = J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1];
    return sqrt(std::abs(symDeterminant(B))) * w2D;
  }

  //==========================================================================
  // STIFFNESS — sparse two-lambda form (used by the SLS attenuation path).
  // Voigt ordering matches makutu, so the assembled term is identical.
  //==========================================================================

  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiBGradPhi(real_t const (&B)[6], FUNC1 &&func1, FUNC2 &&func2) {
    const real_t w = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc));
    func1(qa, qb, qc);
    constexpr int rp1 = num1dNodes;
    for (int i = 0; i < rp1; i++) {
      const int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
      const int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
      const int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
      const real_t gia = basisGradientAt(i, qa);
      const real_t gib = basisGradientAt(i, qb);
      const real_t gic = basisGradientAt(i, qc);
      for (int j = 0; j < rp1; j++) {
        const int jbc = GL_BASIS::TensorProduct3D::linearIndex(j, qb, qc);
        const int ajc = GL_BASIS::TensorProduct3D::linearIndex(qa, j, qc);
        const int abj = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, j);
        const real_t gja = basisGradientAt(j, qa);
        const real_t gjb = basisGradientAt(j, qb);
        const real_t gjc = basisGradientAt(j, qc);
        const real_t w0 = w * gia * gja;
        func2(ibc, jbc, w0 * B[0]);
        const real_t w1 = w * gib * gjb;
        func2(aic, ajc, w1 * B[1]);
        const real_t w2 = w * gic * gjc;
        func2(abi, abj, w2 * B[2]);
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

  /**
   * @brief Sparse stiffness assembly (two-lambda makutu API).
   * @param X     8 corner coordinates.
   * @param func1 func1(qa,qb,qc) — invoked once per quad point (model-on-nodes).
   * @param func2 func2(i,j,R_ij) — invoked per stiffness contribution.
   */
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

  //==========================================================================
  // ACOUSTIC STIFFNESS via sum-factorization — FLAT (RangePolicy) path.
  //
  // Register-lean (3 flux buffers), bit-equivalent to makutu. The physics
  // factor alpha (= 1/rho) is folded into the flux scaling, exactly as in
  // makutu (scale = w * alpha). This is what the default acoustic dispatch
  // (_Acoustic_Flat) calls.
  //==========================================================================

  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(float const (&X)[8][3], real_t const (&u_local)[numNodes],
                                                            real_t (&v_local)[numNodes], FUNC_ALPHA &&get_alpha) {
    real_t G_xi[numNodes] = {0};
    real_t G_eta[numNodes] = {0};
    real_t G_zeta[numNodes] = {0};

    // Pass 1+2: reference gradient of u, then metric + alpha + weight.
    triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;
      constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
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

      real_t const scale = w * get_alpha(qa, qb, qc);  // alpha = 1/rho for acoustics
      G_xi[q] = scale * (B[0] * dxi_q + B[5] * deta_q + B[4] * dzeta_q);
      G_eta[q] = scale * (B[5] * dxi_q + B[1] * deta_q + B[3] * dzeta_q);
      G_zeta[q] = scale * (B[4] * dxi_q + B[3] * deta_q + B[2] * dzeta_q);
    });

    // Pass 3: divergence  v += D^T G^xi + D^T G^eta + D^T G^zeta
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

#ifdef USE_KOKKOS
  //==========================================================================
  // ACOUSTIC STIFFNESS via sum-factorization — TEAM path matching the solver's
  // EXISTING scratch contract (G_xi/G_eta/G_zeta of size numNodes each).
  //
  // One block per element, work distributed over TeamThreadRange. This is the
  // method reached by computeElementContributions_Acoustic_Teams. It fits the
  // solver's current 5*kPtsPerElem scratch budget. alpha is folded into scale.
  //==========================================================================

  template <typename TEAM_MEMBER, typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact_Team(const TEAM_MEMBER &team, float const (&X)[8][3],
                                                                 real_t const *p_local, real_t *f_local, real_t *G_xi,
                                                                 real_t *G_eta, real_t *G_zeta,
                                                                 FUNC_ALPHA &&get_alpha) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numNodes), [&](const int q) {
      int qa, qb, qc;
      GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
      real_t const w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);

      real_t dxi_q = 0, deta_q = 0, dzeta_q = 0;
#pragma unroll
      for (int i = 0; i < num1dNodes; ++i) {
        int const ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
        int const aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
        int const abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
        dxi_q += basisGradientAt(i, qa) * p_local[ibc];
        deta_q += basisGradientAt(i, qb) * p_local[aic];
        dzeta_q += basisGradientAt(i, qc) * p_local[abi];
      }

      real_t J[3][3] = {{0}};
      real_t B[6] = {0};
      computeBMatrix(qa, qb, qc, X, J, B);
      real_t const scale = w * get_alpha(qa, qb, qc);

      G_xi[q] = scale * (B[0] * dxi_q + B[5] * deta_q + B[4] * dzeta_q);
      G_eta[q] = scale * (B[5] * dxi_q + B[1] * deta_q + B[3] * dzeta_q);
      G_zeta[q] = scale * (B[4] * dxi_q + B[3] * deta_q + B[2] * dzeta_q);
    });

    team.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numNodes), [&](const int node) {
      int ia, ib, ic;
      GL_BASIS::TensorProduct3D::multiIndex(node, ia, ib, ic);
      real_t v = 0;
#pragma unroll
      for (int qa = 0; qa < num1dNodes; ++qa) {
        int const q_xi = GL_BASIS::TensorProduct3D::linearIndex(qa, ib, ic);
        v += basisGradientAt(ia, qa) * G_xi[q_xi];
      }
#pragma unroll
      for (int qb = 0; qb < num1dNodes; ++qb) {
        int const q_eta = GL_BASIS::TensorProduct3D::linearIndex(ia, qb, ic);
        v += basisGradientAt(ib, qb) * G_eta[q_eta];
      }
#pragma unroll
      for (int qc = 0; qc < num1dNodes; ++qc) {
        int const q_zeta = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, qc);
        v += basisGradientAt(ic, qc) * G_zeta[q_zeta];
      }
      f_local[node] += v;
    });

    team.team_barrier();
  }

  //==========================================================================
  // ====================  THIRD-PARTY GEMM ACCELERATION PATH  ================
  //
  // Preserved UNCHANGED from the original tensorial implementation. This is
  // where the high-order (> 3) speedup lives: the O(n^4) work is recast as
  // dense [N]x[N^2] GEMMs staged in TEAM SCRATCH (shared) memory, instead of
  // one thread per element holding O(n^3) data in registers.
  //
  // IMPORTANT — these are NOT reached by the API entry points above:
  //   * They require `scratchBytesPerTeam[Streaming]()` bytes of team scratch,
  //     which is LARGER than the solver's current `_Teams` allocation
  //     (5*kPtsPerElem). Wiring them in requires a small solver-side change:
  //     set the policy scratch to scratchBytesPerTeam*() and call
  //     computeStiffnessOperatorTeamVector[Streaming] from a dedicated
  //     computeElementContributions_Acoustic_Gemm kernel.
  //   * The streaming W = w*B does NOT include alpha (1/rho). If you route the
  //     acoustic operator through these, fold 1/rho into W_local (multiply each
  //     6-tuple by get_alpha(qa,qb,qc)) — otherwise the density term is lost.
  //==========================================================================

  /// Serial small matmul: C = A * B  (NN).
  template <int ROWS, int INNER, int COLS>
  PROXY_HOST_DEVICE static void matmul_NN(real_t const (&A)[ROWS][INNER], real_t const (&B)[INNER][COLS],
                                          real_t (&C)[ROWS][COLS]) {
    for (int i = 0; i < ROWS; ++i)
      for (int j = 0; j < COLS; ++j) {
        real_t sum = real_t(0);
        for (int k = 0; k < INNER; ++k) sum += A[i][k] * B[k][j];
        C[i][j] = sum;
      }
  }

  /// Serial small matmul: C = A^T * B  (TN).
  template <int ROWS, int INNER, int COLS>
  PROXY_HOST_DEVICE static void matmul_TN(real_t const (&A)[INNER][ROWS], real_t const (&B)[INNER][COLS],
                                          real_t (&C)[ROWS][COLS]) {
    for (int i = 0; i < ROWS; ++i)
      for (int j = 0; j < COLS; ++j) {
        real_t sum = real_t(0);
        for (int k = 0; k < INNER; ++k) sum += A[k][i] * B[k][j];
        C[i][j] = sum;
      }
  }

  /**
   * @brief Serial (one-thread-per-element) GEMM stiffness operator.
   * @param u  element-local input field (size numNodes)
   * @param Y  element-local output field (size numNodes)  [overwritten]
   * @param W  pre-computed weighted metrics (numNodes*6)
   * @param D  1D reference gradient operator [num1dNodes][num1dNodes],
   *           with D[row][col] == basisGradientAt(col,row).
   */
  PROXY_HOST_DEVICE
  static void computeStiffnessOperatorDevice(real_t const *u, real_t *Y, real_t const *W,
                                             real_t const (&D)[num1dNodes][num1dNodes]) {
    real_t U_Xi[num1dNodes][numNodesPerFace] = {{0}};
    real_t U_Eta[num1dNodes][numNodesPerFace] = {{0}};
    real_t U_Zeta[num1dNodes][numNodesPerFace] = {{0}};

    for (int q_local = 0; q_local < numNodes; q_local++) {
      const int qa = q_local % num1dNodes;
      const int tmp = q_local / num1dNodes;
      const int qb = tmp % num1dNodes;
      const int qc = tmp / num1dNodes;
      U_Xi[qa][qb * num1dNodes + qc] = u[q_local];
      U_Eta[qb][qa * num1dNodes + qc] = u[q_local];
      U_Zeta[qc][qa * num1dNodes + qb] = u[q_local];
    }

    real_t dU_Xi[num1dNodes][numNodesPerFace];
    real_t dU_Eta[num1dNodes][numNodesPerFace];
    real_t dU_Zeta[num1dNodes][numNodesPerFace];
    matmul_NN<num1dNodes, num1dNodes, numNodesPerFace>(D, U_Xi, dU_Xi);
    matmul_NN<num1dNodes, num1dNodes, numNodesPerFace>(D, U_Eta, dU_Eta);
    matmul_NN<num1dNodes, num1dNodes, numNodesPerFace>(D, U_Zeta, dU_Zeta);

    real_t Fx[num1dNodes][numNodesPerFace];
    real_t Fy[num1dNodes][numNodesPerFace];
    real_t Fz[num1dNodes][numNodesPerFace];
    for (int q_local = 0; q_local < numNodes; q_local++) {
      const int qa = q_local % num1dNodes;
      const int tmp = q_local / num1dNodes;
      const int qb = tmp % num1dNodes;
      const int qc = tmp / num1dNodes;
      const int col_Xi = qb * num1dNodes + qc;
      const int col_Eta = qa * num1dNodes + qc;
      const int col_Zeta = qa * num1dNodes + qb;
      const int w_offset = q_local * 6;
      const real_t W0 = W[w_offset + 0], W1 = W[w_offset + 1], W2 = W[w_offset + 2];
      const real_t W3 = W[w_offset + 3], W4 = W[w_offset + 4], W5 = W[w_offset + 5];
      Fx[qa][col_Xi] = W0 * dU_Xi[qa][col_Xi] + W5 * dU_Eta[qb][col_Eta] + W4 * dU_Zeta[qc][col_Zeta];
      Fy[qb][col_Eta] = W5 * dU_Xi[qa][col_Xi] + W1 * dU_Eta[qb][col_Eta] + W3 * dU_Zeta[qc][col_Zeta];
      Fz[qc][col_Zeta] = W4 * dU_Xi[qa][col_Xi] + W3 * dU_Eta[qb][col_Eta] + W2 * dU_Zeta[qc][col_Zeta];
    }

    real_t y_Xi[num1dNodes][numNodesPerFace];
    real_t y_Eta[num1dNodes][numNodesPerFace];
    real_t y_Zeta[num1dNodes][numNodesPerFace];
    matmul_TN<num1dNodes, num1dNodes, numNodesPerFace>(D, Fx, y_Xi);
    matmul_TN<num1dNodes, num1dNodes, numNodesPerFace>(D, Fy, y_Eta);
    matmul_TN<num1dNodes, num1dNodes, numNodesPerFace>(D, Fz, y_Zeta);

    for (int q_local = 0; q_local < numNodes; q_local++) {
      const int qa = q_local % num1dNodes;
      const int tmp = q_local / num1dNodes;
      const int qb = tmp % num1dNodes;
      const int qc = tmp / num1dNodes;
      const int col_Xi = qb * num1dNodes + qc;
      const int col_Eta = qa * num1dNodes + qc;
      const int col_Zeta = qa * num1dNodes + qb;
      Y[q_local] = y_Xi[qa][col_Xi] + y_Eta[qb][col_Eta] + y_Zeta[qc][col_Zeta];
    }
  }

  /// Scratch needed per team for the GEMM team operator (12 [N][N^2] tensors).
  static constexpr size_t scratchBytesPerTeam() {
    constexpr size_t sva = (sizeof(real_t) >= 8) ? sizeof(real_t) : 8;
    constexpr size_t per2d = num1dNodes * numNodesPerFace * sizeof(real_t) + sva;
    return 12 * per2d;
  }

  /// Scratch for the streaming variant (adds a W_local of numNodes*6 reals).
  static constexpr size_t scratchBytesPerTeamStreaming() {
    constexpr size_t sva = (sizeof(real_t) >= 8) ? sizeof(real_t) : 8;
    constexpr size_t perW = numNodes * 6 * sizeof(real_t) + sva;
    return scratchBytesPerTeam() + perW;
  }

  /**
   * @brief Fill the 1D reference-gradient operator expected by the GEMM ops.
   *
   * Row-major, size num1dNodes*num1dNodes, with the convention required by
   * matmul_NN/TN_team: D_flat[row*num1dNodes + col] == basisGradientAt(col,row),
   * i.e. (D u)_row = sum_col phi'_col(xi_row) * u_col. Call once on the host and
   * deep_copy into a device View before launching the GEMM kernel.
   */
  PROXY_HOST_DEVICE
  static void fillDerivativeMatrix(real_t *D_flat) {
    for (int row = 0; row < num1dNodes; ++row)
      for (int col = 0; col < num1dNodes; ++col) D_flat[row * num1dNodes + col] = basisGradientAt(col, row);
  }

  /**
   * @brief Precompute the weighted metric W = w * alpha * B for one element.
   *
   * Writes numNodes*6 reals into @p W_out (Voigt order [xx,yy,zz,yz,xz,xy]).
   * This is the branchy part (computeBMatrix / basisGradientAt / symInvert) that
   * the streaming kernel re-did every timestep; precompute it ONCE per element
   * and the hot stiffness kernel becomes a pure matmul reading W
   * (computeStiffnessOperatorTeamVector). For acoustics, get_alpha = 1/rho.
   *
   * Memory cost of storing W for the whole mesh: nElements * numNodes * 6 reals.
   */
  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeElementMetrics(float const (&X)[8][3], FUNC_ALPHA &&get_alpha, real_t *W_out) {
    for (int q = 0; q < numNodes; ++q) {
      int qa, qb, qc;
      GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
      real_t J[3][3] = {{0}};
      real_t B[6] = {0};
      computeBMatrix(qa, qb, qc, X, J, B);
      const real_t scale = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc) * get_alpha(qa, qb, qc);
      for (int c = 0; c < 6; ++c) W_out[q * 6 + c] = scale * B[c];
    }
  }

  /// Team-parallel C = A * B (NN), output distributed over TeamVectorRange.
  template <int ROWS, int INNER, int COLS, typename MemberType, typename ViewA, typename ViewB, typename ViewC>
  KOKKOS_INLINE_FUNCTION static void matmul_NN_team(const MemberType &member, const ViewA &A, const ViewB &B,
                                                    const ViewC &C) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, ROWS * COLS), [&](int ij) {
      const int row = ij / COLS;
      const int col = ij % COLS;
      real_t sum = real_t(0);
      for (int k = 0; k < INNER; ++k) sum += A(row, k) * B(k, col);
      C(row, col) = sum;
    });
  }

  /// Team-parallel C = A^T * B (TN), output distributed over TeamVectorRange.
  template <int ROWS, int INNER, int COLS, typename MemberType, typename ViewA, typename ViewB, typename ViewC>
  KOKKOS_INLINE_FUNCTION static void matmul_TN_team(const MemberType &member, const ViewA &A, const ViewB &B,
                                                    const ViewC &C) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, ROWS * COLS), [&](int ij) {
      const int row = ij / COLS;
      const int col = ij % COLS;
      real_t sum = real_t(0);
      for (int k = 0; k < INNER; ++k) sum += A(k, row) * B(k, col);
      C(row, col) = sum;
    });
  }

  /**
   * @brief TeamPolicy GEMM stiffness operator (no KokkosBlas dependency).
   *        Needs scratchBytesPerTeam() bytes at level 0.
   * @param W      pre-computed weighted metrics (numNodes*6). For acoustics,
   *               fold 1/rho into W before calling (see note above).
   * @param D_flat 1D derivative operator (num1dNodes*num1dNodes, row-major),
   *               D_flat[row*num1dNodes+col] == basisGradientAt(col,row).
   */
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION static void computeStiffnessOperatorTeamVector(const MemberType &member, real_t *u, real_t *Y,
                                                                        real_t const *W, real_t const *D_flat) {
    constexpr int n = num1dNodes;
    constexpr int n2 = numNodesPerFace;
    constexpr int nTotal = numNodes;

    using ScratchSpace = typename MemberType::execution_space::scratch_memory_space;
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

    using ConstMatView = Kokkos::View<const real_t **, Kokkos::LayoutRight, Kokkos::AnonymousSpace,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    ConstMatView D_v(D_flat, n, n);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n * n2), [&](int idx) {
      int r = idx / n2;
      int c = idx % n2;
      U_Xi(r, c) = real_t(0);
      U_Eta(r, c) = real_t(0);
      U_Zeta(r, c) = real_t(0);
    });
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;
      U_Xi(qa, qb * n + qc) = u[q];
      U_Eta(qb, qa * n + qc) = u[q];
      U_Zeta(qc, qa * n + qb) = u[q];
    });
    member.team_barrier();

    matmul_NN_team<n, n, n2>(member, D_v, U_Xi, dU_Xi);
    member.team_barrier();
    matmul_NN_team<n, n, n2>(member, D_v, U_Eta, dU_Eta);
    member.team_barrier();
    matmul_NN_team<n, n, n2>(member, D_v, U_Zeta, dU_Zeta);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;
      const int cXi = qb * n + qc;
      const int cEta = qa * n + qc;
      const int cZeta = qa * n + qb;
      const int w_off = q * 6;
      const real_t W0 = W[w_off], W1 = W[w_off + 1], W2 = W[w_off + 2];
      const real_t W3 = W[w_off + 3], W4 = W[w_off + 4], W5 = W[w_off + 5];
      Fx(qa, cXi) = W0 * dU_Xi(qa, cXi) + W5 * dU_Eta(qb, cEta) + W4 * dU_Zeta(qc, cZeta);
      Fy(qb, cEta) = W5 * dU_Xi(qa, cXi) + W1 * dU_Eta(qb, cEta) + W3 * dU_Zeta(qc, cZeta);
      Fz(qc, cZeta) = W4 * dU_Xi(qa, cXi) + W3 * dU_Eta(qb, cEta) + W2 * dU_Zeta(qc, cZeta);
    });
    member.team_barrier();

    matmul_TN_team<n, n, n2>(member, D_v, Fx, y_Xi);
    member.team_barrier();
    matmul_TN_team<n, n, n2>(member, D_v, Fy, y_Eta);
    member.team_barrier();
    matmul_TN_team<n, n, n2>(member, D_v, Fz, y_Zeta);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;
      Y[q] = y_Xi(qa, qb * n + qc) + y_Eta(qb, qa * n + qc) + y_Zeta(qc, qa * n + qb);
    });
    member.team_barrier();  // ensure all Y writes are visible before the caller scatters
  }

  /**
   * @brief Streaming GEMM operator: builds W = w*alpha*B in scratch on the fly,
   *        then delegates to computeStiffnessOperatorTeamVector.
   *        Needs scratchBytesPerTeamStreaming() bytes at level 0.
   *
   * @param get_alpha  Callback get_alpha(qa,qb,qc) -> real_t. For acoustics this
   *                   is 1/rho at the quadrature point; the density factor is
   *                   folded into W exactly as in the makutu sum-fact kernel.
   */
  template <typename MemberType, typename FUNC_ALPHA>
  KOKKOS_INLINE_FUNCTION static void computeStiffnessOperatorTeamVectorStreaming(const MemberType &member, real_t *u,
                                                                                 real_t *Y,
                                                                                 real_t const (&cornerCoords)[8][3],
                                                                                 real_t const *D_flat,
                                                                                 FUNC_ALPHA &&get_alpha) {
    constexpr int n = num1dNodes;
    constexpr int nTotal = numNodes;

    using ScratchSpace = typename MemberType::execution_space::scratch_memory_space;
    using ScratchView1D =
        Kokkos::View<real_t *, Kokkos::LayoutRight, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    ScratchView1D W_local(member.team_scratch(0), nTotal * 6);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nTotal), [&](int q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;
      real_t J[3][3] = {{0}};
      real_t B[6] = {0};
      computeBMatrix(qa, qb, qc, cornerCoords, J, B);
      const real_t scale = quadratureWeight(qa) * quadratureWeight(qb) * quadratureWeight(qc) * get_alpha(qa, qb, qc);
      for (int c = 0; c < 6; ++c) W_local(q * 6 + c) = scale * B[c];
    });
    member.team_barrier();

    computeStiffnessOperatorTeamVector(member, u, Y, W_local.data(), D_flat);
  }
#endif  // USE_KOKKOS

  //==========================================================================
  // Virtual compatibility hooks (match makutu).
  //==========================================================================
  PROXY_HOST_DEVICE virtual int getNumQuadraturePoints() { return numQuadraturePoints; }
  PROXY_HOST_DEVICE virtual int getNumSupportPoints() { return numNodes; }
  PROXY_HOST_DEVICE virtual int getMaxSupportPoints() const { return maxSupportPoints; }

  PROXY_HOST_DEVICE
  ~Qk_Hexahedron_Tensorial_GEMM() = default;
};

//============================================================================
// Per-order aliases
//============================================================================
using Q1_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis1>;
using Q2_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis2>;
using Q3_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis3GL>;
using Q4_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis4GL>;
using Q5_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis5GL>;
using Q6_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis6GL>;
using Q7_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis7GL>;
using Q8_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis8GL>;
using Q9_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis9GL>;

//============================================================================
// Order selector (mirrors makutu's selector convention)
//============================================================================
template <int ORDER>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector;

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<1> {
  using type = Q1_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<2> {
  using type = Q2_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<3> {
  using type = Q3_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<4> {
  using type = Q4_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<5> {
  using type = Q5_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<6> {
  using type = Q6_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<7> {
  using type = Q7_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<8> {
  using type = Q8_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM_Selector<9> {
  using type = Q9_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM;
};

#endif  // FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_TENSORIAL_GEMM_H_
