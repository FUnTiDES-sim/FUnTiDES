#ifndef SRC_DISCRETIZATION_FE_MAKUTU_INCLUDE_QK_HEXAHEDRON_TENSOSIRAL_H_
#define SRC_DISCRETIZATION_FE_MAKUTU_INCLUDE_QK_HEXAHEDRON_TENSOSIRAL_H_

#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"

template <typename GL_BASIS>
class Qk_Hexahedron_Tensorial : public Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS> {
 public:
  using Base = Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>;

  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;
  constexpr static int numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;
  constexpr static int numNodesPerFace = num1dNodes * num1dNodes;

  template <int ROWS, int INNER, int COLS>
  PROXY_HOST_DEVICE static void matmul_NN(real_t const (&A)[ROWS][INNER], real_t const (&B)[INNER][COLS],
                                          real_t (&C)[ROWS][COLS]) {
    for (int i = 0; i < ROWS; ++i) {
      for (int j = 0; j < COLS; ++j) {
        real_t sum = real_t(0);
        for (int k = 0; k < INNER; ++k) sum += A[i][k] * B[k][j];
        C[i][j] = sum;
      }
    }
  }

  template <int ROWS, int INNER, int COLS>
  PROXY_HOST_DEVICE static void matmul_TN(real_t const (&A)[INNER][ROWS], real_t const (&B)[INNER][COLS],
                                          real_t (&C)[ROWS][COLS]) {
    for (int i = 0; i < ROWS; ++i) {
      for (int j = 0; j < COLS; ++j) {
        real_t sum = real_t(0);
        for (int k = 0; k < INNER; ++k) sum += A[k][i] * B[k][j];
        C[i][j] = sum;
      }
    }
  }

  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(float const (&X)[8][3], real_t const (&p_local)[numNodes],
                                                            real_t (&f_local)[numNodes], FUNC_ALPHA &&get_alpha) {
    constexpr int n = num1dNodes;
    constexpr int n2 = numNodesPerFace;

    real_t D[n][n];
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        D[i][j] = Base::basisGradientAt(j, i);
      }
    }

    real_t U_Xi[n][n2] = {{0}}, U_Eta[n][n2] = {{0}}, U_Zeta[n][n2] = {{0}};
    for (int q = 0; q < numNodes; ++q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;
      U_Xi[qa][qb * n + qc] = p_local[q];
      U_Eta[qb][qa * n + qc] = p_local[q];
      U_Zeta[qc][qa * n + qb] = p_local[q];
    }

    real_t dU_Xi[n][n2], dU_Eta[n][n2], dU_Zeta[n][n2];
    matmul_NN<n, n, n2>(D, U_Xi, dU_Xi);
    matmul_NN<n, n, n2>(D, U_Eta, dU_Eta);
    matmul_NN<n, n, n2>(D, U_Zeta, dU_Zeta);

    real_t Fx[n][n2], Fy[n][n2], Fz[n][n2];
    for (int q = 0; q < numNodes; ++q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;

      real_t J[3][3] = {{0}};
      real_t B[6] = {0};
      Base::computeBMatrix(qa, qb, qc, X, J, B);

      const real_t w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);
      const real_t scale = w * get_alpha(qa, qb, qc);

      const real_t W0 = scale * B[0];
      const real_t W1 = scale * B[1];
      const real_t W2 = scale * B[2];
      const real_t W3 = scale * B[3];
      const real_t W4 = scale * B[4];
      const real_t W5 = scale * B[5];

      const int cXi = qb * n + qc;
      const int cEta = qa * n + qc;
      const int cZeta = qa * n + qb;

      Fx[qa][cXi] = W0 * dU_Xi[qa][cXi] + W5 * dU_Eta[qb][cEta] + W4 * dU_Zeta[qc][cZeta];
      Fy[qb][cEta] = W5 * dU_Xi[qa][cXi] + W1 * dU_Eta[qb][cEta] + W3 * dU_Zeta[qc][cZeta];
      Fz[qc][cZeta] = W4 * dU_Xi[qa][cXi] + W3 * dU_Eta[qb][cEta] + W2 * dU_Zeta[qc][cZeta];
    }

    real_t y_Xi[n][n2], y_Eta[n][n2], y_Zeta[n][n2];
    matmul_TN<n, n, n2>(D, Fx, y_Xi);
    matmul_TN<n, n, n2>(D, Fy, y_Eta);
    matmul_TN<n, n, n2>(D, Fz, y_Zeta);

    for (int q = 0; q < numNodes; ++q) {
      const int qa = q % n;
      const int tmp = q / n;
      const int qb = tmp % n;
      const int qc = tmp / n;
      f_local[q] += y_Xi[qa][qb * n + qc] + y_Eta[qb][qa * n + qc] + y_Zeta[qc][qa * n + qb];
    }
  }
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
