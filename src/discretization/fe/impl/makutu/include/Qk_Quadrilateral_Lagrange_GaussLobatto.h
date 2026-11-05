/**
 * @file Qk_Quadrilateral_Lagrange_GaussLobatto.h
 * @brief Element kernels of the makutu back-end in 2-D: Qk quadrilaterals with
 * Gauss-Lobatto-Legendre nodes, and the compile-time loop helper double_loop().
 */

#ifndef FUNTIDES_DISCRETIZATION_FE_IMPL_MAKUTU_INCLUDE_QK_QUADRILATERAL_LAGRANGE_GAUSSLOBATTO_H_
#define FUNTIDES_DISCRETIZATION_FE_IMPL_MAKUTU_INCLUDE_QK_QUADRILATERAL_LAGRANGE_GAUSSLOBATTO_H_

#include <data_type.h>

#include <cmath>

// 1D Lagrange bases (LagrangeBasis1 ... LagrangeBasis9GL) and for_constexpr()
#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"

/**
 * @brief Unrolled double loop: calls lambda(I, J) with std::integral_constant
 * arguments for I in [0, BoundI) and J in [0, BoundJ), J varying fastest.
 */
template <int BoundI, int BoundJ, typename Lambda>
constexpr void double_loop(Lambda &&lambda) {
  for_constexpr<BoundI>([&](auto I) { for_constexpr<BoundJ>([&](auto J) { lambda(I, J); }); });
}

/**
 * @brief Spectral-element kernels of a Qk quadrilateral whose nodes, and
 * quadrature points, are the tensor-product Gauss-Lobatto-Legendre points.
 *
 * Two-dimensional counterpart of Qk_Hexahedron_Lagrange_GaussLobatto, for the
 * 2-D solvers in the (x, z) plane. A stateless collection of static host/device
 * functions called one element at a time: local index maps, Jacobian of the
 * bilinear map from the parent square [-1, 1]^2, and the element mass,
 * stiffness and boundary (damping) terms. The mass matrix is diagonal because
 * nodes and quadrature points coincide.
 *
 * Local numbering: node (qa, qb) has index qa + qb * num1dNodes, qa along the
 * first parent axis (x) and qb along the second one (z). The 4 vertices are
 * numbered as the nodes, k % 2 along x and k / 2 along z: X[k][0] = x,
 * X[k][1] = z.
 *
 * Select a degree with Qk_Quadrilateral_Lagrange_GaussLobatto_Selector.
 * @tparam GL_BASIS 1D Lagrange basis, as for the hexahedron (its
 * TensorProduct2D gives the node numbering).
 */
template <typename GL_BASIS>
class Qk_Quadrilateral_Lagrange_GaussLobatto final {
 public:
  /// Space dimension.
  constexpr static int numDims = 2;
  /// Number of nodes along one parent axis (order + 1).
  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;
  /// Half of (num1dNodes - 1), used by basisGradientAt().
  constexpr static int halfNodes = (GL_BASIS::numSupportPoints - 1) / 2;
  /// Number of nodes of the element, num1dNodes^2.
  constexpr static int numNodes = GL_BASIS::TensorProduct2D::numSupportPoints;
  /// Number of nodes on one edge.
  constexpr static int numNodesPerEdge = num1dNodes;
  /// Number of quadrature points (GLL collocation: one per node).
  constexpr static int numQuadraturePoints = numNodes;
  /// Number of vertices.
  constexpr static int numVertices = 4;

  /**
   * @brief Element-local index of node (qa, qb).
   */
  PROXY_HOST_DEVICE
  constexpr static int linearIndex2D(int const qa, int const qb) { return GL_BASIS::TensorProduct2D::linearIndex(qa, qb); }

  /**
   * @brief Element-local index of vertex k (0 to 3).
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex2D(int const k) {
    return linearIndex2D((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * (k / 2));
  }

  /**
   * @brief Element-local index of node i (0 to num1dNodes - 1) of an edge.
   *
   * Edges: 0 is z = zmin (qb = 0), 1 is x = xmax (qa = num1dNodes - 1), 2 is
   * z = zmax (qb = num1dNodes - 1), 3 is x = xmin (qa = 0). Node i runs along
   * the edge in the increasing direction of the other parent axis.
   */
  PROXY_HOST_DEVICE
  constexpr static int edgeNodeIndex(int const edge, int const i) {
    constexpr int last = num1dNodes - 1;
    return edge == 0 ? linearIndex2D(i, 0) : edge == 1 ? linearIndex2D(last, i) : edge == 2 ? linearIndex2D(i, last) : linearIndex2D(0, i);
  }

  /**
   * @brief Parent coordinate of the 1D node q, in [-1, 1].
   */
  PROXY_HOST_DEVICE
  constexpr static real_t parentCoord(int const q) { return static_cast<real_t>(GL_BASIS::parentSupportCoord(q)); }

  /**
   * @brief 1D GLL quadrature weight of node q.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t quadratureWeight(int const q) { return static_cast<real_t>(GL_BASIS::weight(q)); }

  /**
   * @brief Linear interpolation coefficient of the interval end k (0 or 1) at
   * node q: (1 - xi_q) / 2 for k = 0, (1 + xi_q) / 2 for k = 1.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t interpolationCoord(int const q, int const k) {
    const real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0);
    return k == 0 ? (1.0 - alpha) : alpha;
  }

  /**
   * @brief Derivative of the 1D basis function q at node p.
   *
   * Uses the antisymmetry of the GLL derivative table about the center.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t basisGradientAt(int const q, int const p) {
    if (p <= halfNodes) {
      return static_cast<real_t>(GL_BASIS::gradientAt(q, p));
    }
    return -static_cast<real_t>(GL_BASIS::gradientAt(num1dNodes - 1 - q, num1dNodes - 1 - p));
  }

  /**
   * @brief Shape function values at parent coordinates.
   * @param[in] coords Parent coordinates (xi, eta) in [-1, 1]^2.
   * @param[out] N Value of each element node's shape function.
   */
  PROXY_HOST_DEVICE
  static void calcN(double const (&coords)[2], double (&N)[numNodes]) { GL_BASIS::TensorProduct2D::value(coords, N); }

  /**
   * @brief Jacobian of the bilinear map at node (qa, qb): J[i][j] = d x_i / d xi_j.
   * @param[in] X Coordinates of the 4 vertices, X[vertex][axis].
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, real_t const (&X)[4][2], real_t (&J)[2][2]) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j) J[i][j] = 0;
    for (int k = 0; k < 4; ++k) {
      int const ka = k % 2;
      int const kb = k / 2;
      // d/dxi of the vertex shape function, then d/deta
      real_t const dxi = (ka == 0 ? -0.5 : 0.5) * interpolationCoord(qb, kb);
      real_t const deta = interpolationCoord(qa, ka) * (kb == 0 ? -0.5 : 0.5);
      for (int i = 0; i < 2; ++i) {
        J[i][0] += X[k][i] * dxi;
        J[i][1] += X[k][i] * deta;
      }
    }
  }

  /**
   * @brief Determinant of a 2x2 matrix.
   */
  PROXY_HOST_DEVICE
  static real_t determinant2x2(real_t const (&J)[2][2]) { return J[0][0] * J[1][1] - J[0][1] * J[1][0]; }

  /**
   * @brief Stiffness metric B = |det J| J^-1 J^-T at node (qa, qb).
   * @param[in] X Coordinates of the 4 vertices.
   * @param[out] J Jacobian at the node.
   * @param[out] B Voigt order [xi xi, eta eta, xi eta].
   * @return det(J).
   */
  PROXY_HOST_DEVICE
  static real_t computeBMatrix(int const qa, int const qb, real_t const (&X)[4][2], real_t (&J)[2][2], real_t (&B)[3]) {
    jacobianTransformation(qa, qb, X, J);
    real_t const detJ = determinant2x2(J);
    // J^-1 = adj(J) / det(J), so |det J| J^-1 J^-T = adj(J) adj(J)^T / |det J|
    real_t const inv = 1.0 / std::abs(detJ);
    B[0] = (J[1][1] * J[1][1] + J[0][1] * J[0][1]) * inv;
    B[1] = (J[1][0] * J[1][0] + J[0][0] * J[0][0]) * inv;
    B[2] = -(J[1][1] * J[1][0] + J[0][1] * J[0][0]) * inv;
    return detJ;
  }

  /**
   * @brief Diagonal of the element mass matrix for a unit coefficient:
   * |det J| times the quadrature weight, at each node.
   * @param[in] X Coordinates of the 4 vertices.
   * @param[in] func Called as func(q, value) once per node.
   */
  template <typename FUNC>
  PROXY_HOST_DEVICE static void computeMassTerm(real_t const (&X)[4][2], FUNC &&func) {
    double_loop<num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int q = GL_BASIS::TensorProduct2D::linearIndex(qa, qb);
      constexpr real_t w2D = GL_BASIS::weight(qa) * GL_BASIS::weight(qb);
      real_t J[2][2];
      jacobianTransformation(qa, qb, X, J);
      func(q, std::abs(determinant2x2(J)) * w2D);
    });
  }

  /**
   * @brief Adds the element stiffness operator applied to @p u_local to
   * @p v_local, by sum factorization (O(order^3) per element).
   *
   * v_local[i] += sum_q w_q alpha_q (grad phi_i)^T B_q (grad u)(q), with parent
   * gradients and B from computeBMatrix().
   * @param[in] X Coordinates of the 4 vertices.
   * @param[in] u_local Nodal values, element-local numbering.
   * @param[in,out] v_local Result; must not alias @p u_local.
   * @param[in] get_alpha Called as get_alpha(qa, qb); coefficient at the node
   * (e.g. 1 / rho for the acoustic equation).
   */
  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(real_t const (&X)[4][2], real_t const (&u_local)[numNodes],
                                                            real_t (&v_local)[numNodes], FUNC_ALPHA &&get_alpha) {
    // Weighted parent-coordinate fluxes at each quadrature point
    real_t G_xi[numNodes] = {0};
    real_t G_eta[numNodes] = {0};

    double_loop<num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int q = GL_BASIS::TensorProduct2D::linearIndex(qa, qb);
      constexpr real_t w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb);

      real_t dxi = 0, deta = 0;
      for_constexpr<num1dNodes>([&](auto ici) {
        constexpr int i = decltype(ici)::value;
        dxi += basisGradientAt(i, qa) * u_local[GL_BASIS::TensorProduct2D::linearIndex(i, qb)];
        deta += basisGradientAt(i, qb) * u_local[GL_BASIS::TensorProduct2D::linearIndex(qa, i)];
      });

      real_t J[2][2];
      real_t B[3];
      computeBMatrix(qa, qb, X, J, B);
      real_t const scale = w * get_alpha(qa, qb);
      G_xi[q] = scale * (B[0] * dxi + B[2] * deta);
      G_eta[q] = scale * (B[2] * dxi + B[1] * deta);
    });

    // v[node] += sum over the quadrature points on the two lines through the node
    double_loop<num1dNodes, num1dNodes>([&](auto const icia, auto const icib) {
      constexpr int ia = decltype(icia)::value;
      constexpr int ib = decltype(icib)::value;
      real_t v = 0;
      for_constexpr<num1dNodes>([&](auto icq) {
        constexpr int qq = decltype(icq)::value;
        v += basisGradientAt(ia, qq) * G_xi[GL_BASIS::TensorProduct2D::linearIndex(qq, ib)];
        v += basisGradientAt(ib, qq) * G_eta[GL_BASIS::TensorProduct2D::linearIndex(ia, qq)];
      });
      v_local[GL_BASIS::TensorProduct2D::linearIndex(ia, ib)] += v;
    });
  }

  /**
   * @brief Edge measure times quadrature weight at an edge node: the diagonal
   * entry of the edge mass matrix (absorbing-boundary term).
   * @param[in] i Edge-local node index (0 to num1dNodes - 1).
   * @param[in] X Coordinates of the 2 edge vertices, X[vertex][axis]; the edge
   * is straight.
   */
  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const i, real_t const (&X)[2][2]) {
    real_t const dx = X[1][0] - X[0][0];
    real_t const dz = X[1][1] - X[0][1];
    return static_cast<real_t>(GL_BASIS::weight(i)) * 0.5 * std::sqrt(dx * dx + dz * dz);
  }
};

/// @name Quadrilaterals of order 1 to 9
/// @{
using Q1_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis1>;
using Q2_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis2>;
using Q3_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis3GL>;
using Q4_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis4GL>;
using Q5_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis5GL>;
using Q6_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis6GL>;
using Q7_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis7GL>;
using Q8_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis8GL>;
using Q9_Quadrilateral_Lagrange_GaussLobatto = Qk_Quadrilateral_Lagrange_GaussLobatto<LagrangeBasis9GL>;
/// @}

/**
 * @brief Maps a polynomial order (1 to 9) to its quadrilateral, exposed as `type`.
 */
template <int ORDER>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector;

/// @cond Doxygen_Suppress
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<1> {
  using type = Q1_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<2> {
  using type = Q2_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<3> {
  using type = Q3_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<4> {
  using type = Q4_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<5> {
  using type = Q5_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<6> {
  using type = Q6_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<7> {
  using type = Q7_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<8> {
  using type = Q8_Quadrilateral_Lagrange_GaussLobatto;
};
template <>
struct Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<9> {
  using type = Q9_Quadrilateral_Lagrange_GaussLobatto;
};
/// @endcond

#endif  // FUNTIDES_DISCRETIZATION_FE_IMPL_MAKUTU_INCLUDE_QK_QUADRILATERAL_LAGRANGE_GAUSSLOBATTO_H_
