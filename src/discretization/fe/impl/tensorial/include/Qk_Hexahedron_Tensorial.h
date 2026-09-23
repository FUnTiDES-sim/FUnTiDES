#ifndef FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_TENSORIAL_GEMM_H_
#define FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_TENSORIAL_GEMM_H_

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
#include "mathUtilites.h"

/**
 * @brief Element kernels of a Qk Gauss-Lobatto hexahedron (mass, damping, stiffness).
 *
 * Nodes and quadrature points coincide. Local node index is
 * qa + qb * num1dNodes + qc * num1dNodes^2 (qa runs fastest), see linearIndex3DVal.
 * Vertex k of an element (8 corners) has 3D indices (k % 2, (k % 4) / 2, k / 4).
 * Voigt order for symmetric 3x3 tensors is [xx, yy, zz, yz, xz, xy].
 *
 * The stiffness operators are available as a sparse two-callback form, a flat
 * sum-factorization form and dense-GEMM forms (serial and Kokkos team-parallel).
 *
 * @tparam GL_BASIS Lagrange basis type (LagrangeBasis1, LagrangeBasis2, LagrangeBasis3GL, ...).
 */
template <typename GL_BASIS>
class Qk_Hexahedron_Tensorial_GEMM final {
 public:
  using BasisType = GL_BASIS;  ///< Basis type, exposed for tests and external use.

  /// Number of nodes per element per dimension.
  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;

  /// Half the number of nodes per dimension, rounded down ((num1dNodes - 1) / 2).
  constexpr static int halfNodes = (GL_BASIS::numSupportPoints - 1) / 2;

  /// Total number of nodes per element.
  constexpr static int numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;

  /// Number of nodes per element face.
  constexpr static int numNodesPerFace = num1dNodes * num1dNodes;

  /// Maximum number of support points per element.
  constexpr static int maxSupportPoints = numNodes;

  /// Number of quadrature points per element.
  constexpr static int numQuadraturePoints = numNodes;

  /// Jacobian storage type.
  struct JacobianType {
    float data[3][3];
  };

  /// Tag type.
  struct TeamGemm {};

  /// @brief Local node index of 3D indices (qa, qb, qc); qa runs fastest.
  PROXY_HOST_DEVICE
  constexpr static int linearIndex3DVal(const int qa, int const qb, int const qc) {
    return qa + qb * num1dNodes + qc * numNodesPerFace;
  }

  /// @brief Local node index of element vertex k (0 to 7), see the class comment for the vertex order.
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex3D(int const k) {
    return linearIndex3DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * ((k % 4) / 2), (num1dNodes - 1) * (k / 4));
  }

  /// @brief Local face node index of 2D indices (qa, qb); qa runs fastest.
  PROXY_HOST_DEVICE
  constexpr static int linearIndex2DVal(const int qa, const int qb) { return qa + qb * num1dNodes; }

  /// @brief Local face node index of face vertex k (0 to 3), with 2D indices (k % 2, k / 2).
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex2D(int const k) {
    return linearIndex2DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * (k / 2));
  }

  /**
   * @brief Derivative of the 1D Lagrange polynomial of node p, evaluated at node q, in reference coordinates.
   * @param q Evaluation node, in [0, num1dNodes).
   * @param p Polynomial node, in [0, num1dNodes).
   * @return d(phi_p)/d(xi) at xi_q.
   *
   * Only the first half of the basis derivative table is used; the second half is
   * obtained by the Gauss-Lobatto symmetry.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t basisGradientAt(const int q, const int p) {
    if (p <= halfNodes) {
      return GL_BASIS::gradientAt(q, p);
    } else {
      return -GL_BASIS::gradientAt(num1dNodes - 1 - q, num1dNodes - 1 - p);
    }
  }

  /// @brief 1D Gauss-Lobatto quadrature weight at node q.
  PROXY_HOST_DEVICE
  constexpr static real_t quadratureWeight(const int q) { return GL_BASIS::weight(q); }

  /**
   * @brief Value at node q of the 1D linear shape function of vertex k.
   * @param q Node index, in [0, num1dNodes).
   * @param k 0 for the vertex at xi = -1, 1 for the vertex at xi = +1.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t interpolationCoord(const int q, const int k) {
    const real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0);
    return k == 0 ? (1.0 - alpha) : alpha;
  }

  /**
   * @brief Factor of the trilinear Jacobian coming from one vertex along one direction.
   * @param q   Node index along direction i.
   * @param i   Direction of this factor (0, 1 or 2).
   * @param k   Vertex index along direction i (0 or 1).
   * @param dir Derivative direction (0, 1 or 2).
   * @return Derivative of the linear shape function if i == dir, its value at node q otherwise.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t jacobianCoefficient1D(const int q, const int i, const int k, const int dir) {
    if (i == dir)
      return k == 0 ? -0.5 : 0.5;
    else
      return interpolationCoord(q, k);
  }

  /**
   * @brief Jacobian of the trilinear (8-vertex) geometric map at node (qa, qb, qc).
   * @param[in]  X 8 vertex coordinates, X[k][i] is coordinate i of vertex k.
   * @param[out] J J[i][j] = d x_i / d xi_j.
   */
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

  /**
   * @brief Jacobian of the bilinear (4-vertex) face geometric map at face node (qa, qb).
   * @param[in]  X 4 face vertex coordinates, X[k][i] is coordinate i of vertex k.
   * @param[out] J J[i][j] = d x_i / d xi_j, j in {0, 1}.
   */
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
   * @brief Metric tensor B = det(J) * (J^T J)^{-1} at node (qa, qb, qc), in Voigt order.
   * @param[in]  X 8 vertex coordinates.
   * @param[out] J Jacobian at the node.
   * @param[out] B Symmetric metric tensor, Voigt order [xx, yy, zz, yz, xz, xy].
   */
  PROXY_HOST_DEVICE
  static void computeBMatrix(int const qa, int const qb, int const qc, real_t const (&X)[8][3], real_t (&J)[3][3],
                             real_t (&B)[6]) {
    jacobianTransformation(qa, qb, qc, X, J);
    real_t const detJ = determinant(J);
    real_t const invDetJ = 1.0 / detJ;

    // B holds J^T J / det(J) here; the inversion below turns it into det(J) * (J^T J)^{-1}.
    B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) * invDetJ;
    B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) * invDetJ;
    B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]) * invDetJ;
    B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]) * invDetJ;
    B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]) * invDetJ;
    B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) * invDetJ;

    symInvert(B);
  }

  /**
   * @brief Diagonal mass contribution |det(J)| * w3D of every node of the element.
   * @tparam FUNC Callable with signature void(int q, real_t val).
   * @param[in] X    8 vertex coordinates.
   * @param     func Called once per node q with its mass contribution val.
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

  /**
   * @brief Diagonal face-integrated (damping) contribution of one face node.
   * @param q Face node index, see linearIndex2DVal.
   * @param X 4 face vertex coordinates.
   * @return Surface factor sqrt(|det(J^T J)|) * w2D.
   */
  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]) {
    int qa, qb;
    GL_BASIS::TensorProduct2D::multiIndex(q, qa, qb);
    const real_t w2D = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb));
    real_t B[3];
    real_t J[3][2] = {{0}};
    jacobianTransformation2d(qa, qb, X, J);
    // J^T J, 2x2 symmetric, Voigt order [00, 11, 01].
    B[0] = J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0];
    B[1] = J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1];
    B[2] = J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1];
    return sqrt(std::abs(symDeterminant(B))) * w2D;
  }

  /**
   * @brief Emits the stiffness contributions of quadrature point (qa, qb, qc).
   * @tparam qa,qb,qc Quadrature point indices.
   * @tparam FUNC1    Callable with signature void(int qa, int qb, int qc).
   * @tparam FUNC2    Callable with signature void(int i, int j, real_t value).
   * @param[in] B     Metric tensor at the point, Voigt order [xx, yy, zz, yz, xz, xy].
   * @param func1     Called once, before any func2 call.
   * @param func2     Called once per contribution to the stiffness entry (i, j), i and j being local node indices.
   */
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
   * @brief Sparse stiffness assembly of one element through two callbacks.
   * @tparam FUNC1 Callable with signature void(int qa, int qb, int qc).
   * @tparam FUNC2 Callable with signature void(int i, int j, real_t value).
   * @param[in] X 8 vertex coordinates.
   * @param func1 Called once per quadrature point, before the func2 calls of that point (used to fetch model values).
   * @param func2 Called once per contribution to the stiffness entry (i, j), i and j being local node indices.
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

  /**
   * @brief Acoustic stiffness of one element by sum-factorization, flat (one thread per element) path.
   *
   * Accumulates the product of the element stiffness matrix and u_local into
   * v_local. The factor alpha returned by get_alpha is folded into the flux scaling.
   *
   * @tparam FUNC_ALPHA Callable with signature real_t(int qa, int qb, int qc).
   * @param[in]     X         8 vertex coordinates.
   * @param[in]     u_local   Nodal input field, size numNodes.
   * @param[in,out] v_local   Nodal output field, size numNodes; the result is added to it.
   * @param         get_alpha Returns the coefficient alpha at quadrature point (qa, qb, qc); 1/rho for acoustics.
   */
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

      real_t const scale = w * get_alpha(qa, qb, qc);
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

  /// @brief Serial small matrix product C = A * B.
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

  /// @brief Serial small matrix product C = A^T * B.
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
   * @brief Stiffness operator of one element by dense GEMMs, executed by a single thread.
   * @param[in]  u Nodal input field, size numNodes.
   * @param[out] Y Nodal output field, size numNodes; overwritten.
   * @param[in]  W Weighted metrics from computeElementMetrics, size numNodes * 6, Voigt order per node.
   * @param[in]  D 1D reference gradient operator, D[row][col] == basisGradientAt(col, row), see fillDerivativeMatrix.
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

  /// @brief Team scratch size in bytes (level 0) needed by computeStiffnessOperatorTeamVector: 12 tensors [num1dNodes][num1dNodes^2].
  static constexpr size_t scratchBytesPerTeam() {
    constexpr size_t sva = (sizeof(real_t) >= 8) ? sizeof(real_t) : 8;
    constexpr size_t per2d = num1dNodes * numNodesPerFace * sizeof(real_t) + sva;
    return 12 * per2d;
  }

  /// @brief Team scratch size in bytes (level 0) needed by computeStiffnessOperatorTeamVectorStreaming: scratchBytesPerTeam() plus numNodes * 6 reals.
  static constexpr size_t scratchBytesPerTeamStreaming() {
    constexpr size_t sva = (sizeof(real_t) >= 8) ? sizeof(real_t) : 8;
    constexpr size_t perW = numNodes * 6 * sizeof(real_t) + sva;
    return scratchBytesPerTeam() + perW;
  }

  /**
   * @brief Fills the 1D reference gradient operator expected by the GEMM operators.
   * @param[out] D_flat Row-major array of size num1dNodes * num1dNodes,
   *                    D_flat[row * num1dNodes + col] == basisGradientAt(col, row).
   *
   * Call once on the host, then copy to a device View before launching a GEMM kernel.
   */
  PROXY_HOST_DEVICE
  static void fillDerivativeMatrix(real_t *D_flat) {
    for (int row = 0; row < num1dNodes; ++row)
      for (int col = 0; col < num1dNodes; ++col) D_flat[row * num1dNodes + col] = basisGradientAt(col, row);
  }

  /**
   * @brief Precomputes the weighted metric W = w3D * alpha * B of one element, for the GEMM operators.
   *
   * Computing W once per element leaves the stiffness kernel as pure matrix
   * products. Storing W for a whole mesh costs nElements * numNodes * 6 reals.
   *
   * @tparam FUNC_ALPHA Callable with signature real_t(int qa, int qb, int qc).
   * @param[in]  X         8 vertex coordinates.
   * @param      get_alpha Coefficient alpha at quadrature point (qa, qb, qc); 1/rho for acoustics.
   * @param[out] W_out     Size numNodes * 6; node q occupies W_out[6*q .. 6*q+5], Voigt order [xx, yy, zz, yz, xz, xy].
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

  /// @brief Team-parallel matrix product C = A * B, output entries distributed over TeamVectorRange. The caller must synchronize the team afterwards.
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

  /// @brief Team-parallel matrix product C = A^T * B, output entries distributed over TeamVectorRange. The caller must synchronize the team afterwards.
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
   * @brief Stiffness operator of one element by dense GEMMs, executed by one Kokkos team.
   *
   * Needs scratchBytesPerTeam() bytes of level-0 team scratch. The output is
   * visible to the team on return.
   *
   * @tparam MemberType Kokkos team member type.
   * @param      member Team handle.
   * @param[in]  u      Nodal input field, size numNodes.
   * @param[out] Y      Nodal output field, size numNodes; overwritten.
   * @param[in]  W      Weighted metrics from computeElementMetrics, size numNodes * 6. Any coefficient (1/rho for acoustics) must already be folded into W.
   * @param[in]  D_flat 1D derivative operator, row-major, see fillDerivativeMatrix.
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
    member.team_barrier();  // all Y writes must be visible before the caller scatters
  }

  /**
   * @brief Team GEMM stiffness operator that builds the weighted metric W in team scratch, then applies it.
   *
   * Needs scratchBytesPerTeamStreaming() bytes of level-0 team scratch. Same
   * contract as computeStiffnessOperatorTeamVector for u, Y and D_flat.
   *
   * @tparam MemberType Kokkos team member type.
   * @tparam FUNC_ALPHA Callable with signature real_t(int qa, int qb, int qc).
   * @param      member       Team handle.
   * @param[in]  u            Nodal input field, size numNodes.
   * @param[out] Y            Nodal output field, size numNodes; overwritten.
   * @param[in]  cornerCoords 8 vertex coordinates.
   * @param[in]  D_flat       1D derivative operator, row-major, see fillDerivativeMatrix.
   * @param      get_alpha    Coefficient alpha at quadrature point (qa, qb, qc); 1/rho for acoustics.
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

  PROXY_HOST_DEVICE virtual int getNumQuadraturePoints() { return numQuadraturePoints; }
  PROXY_HOST_DEVICE virtual int getNumSupportPoints() { return numNodes; }
  PROXY_HOST_DEVICE virtual int getMaxSupportPoints() const { return maxSupportPoints; }

  PROXY_HOST_DEVICE
  ~Qk_Hexahedron_Tensorial_GEMM() = default;
};

using Q1_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis1>;
using Q2_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis2>;
using Q3_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis3GL>;
using Q4_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis4GL>;
using Q5_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis5GL>;
using Q6_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis6GL>;
using Q7_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis7GL>;
using Q8_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis8GL>;
using Q9_Hexahedron_Lagrange_GaussLobatto_Tensorial_GEMM = Qk_Hexahedron_Tensorial_GEMM<LagrangeBasis9GL>;

/**
 * @brief Maps a polynomial order to the corresponding Qk_Hexahedron_Tensorial_GEMM type.
 * @tparam ORDER Polynomial order, 1 to 9. Other values are not defined.
 *
 * The selected type is the member alias `type`.
 */
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
