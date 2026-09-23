#ifndef FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_BASE_H_
#define FUNTIDES_DISCRETIZATION_FE_QK_HEXAHEDRON_BASE_H_

#include <data_type.h>

#include "fe_discretization.h"
#include "mathUtilites.h"

/**
 * @brief Back-end independent part of the Qk hexahedral spectral-element
 * discretization on Gauss-Lobatto nodes.
 *
 * Provides the compile-time sizes, the local index maps, the Jacobian of the
 * trilinear map from the parent cube [-1, 1]^3 to a hexahedron given by its 8
 * vertices, and the element mass, damping and stiffness kernels. Every member
 * is static and host/device callable, since virtual dispatch is unusable on
 * GPU.
 * @tparam GL_BASIS 1D Lagrange basis, with the interface described in
 * docs/design.md, "1D Lagrange bases".
 * @note No back-end derives from this class yet, and it uses triple_loop and
 * for_constexpr without including their definition (see
 * docs/design-red-flags.md).
 * @see docs/design.md, "Hexahedron local numbering".
 */
template <typename GL_BASIS>
class QkHexahedronBase : public discretization::fe::api::FeDiscretizationTag {
 public:
  using BasisType = GL_BASIS;  ///< The 1D Lagrange basis.

  /// Number of nodes per element edge (order + 1).
  constexpr static int num1dNodes = GL_BASIS::numSupportPoints;

  /// Largest 1D node index for which GL_BASIS::gradientAt() is queried; the
  /// other half follows by symmetry.
  constexpr static int halfNodes = (GL_BASIS::numSupportPoints - 1) / 2;

  /// Number of nodes per element, num1dNodes^3.
  constexpr static int numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;

  /// Number of nodes per face, num1dNodes^2.
  constexpr static int numNodesPerFace = GL_BASIS::TensorProduct2D::numSupportPoints;

  /// Number of support points per element, equal to numNodes.
  constexpr static int maxSupportPoints = numNodes;

  /// Number of quadrature points per element: the quadrature points are the
  /// nodes.
  constexpr static int numQuadraturePoints = numNodes;

  /// A 3x3 Jacobian matrix. Unused by this class.
  struct JacobianType {
    float data[3][3];  ///< Matrix entries.
  };

  /**
   * @brief Element-local index of the node (qa, qb, qc).
   * @return qa + qb*num1dNodes + qc*num1dNodes^2.
   * @see docs/design.md, "Hexahedron local numbering".
   */
  PROXY_HOST_DEVICE
  constexpr static int linearIndex3DVal(const int qa, int const qb, int const qc) {
    return qa + qb * num1dNodes + qc * numNodesPerFace;
  }

  /**
   * @brief Element-local index of the node at a vertex of the hexahedron.
   * @param[in] k Vertex index in [0, 7]; bits 0, 1 and 2 give the side (0 =
   * minus, 1 = plus) along the first, second and third parent axes.
   * @return The element-local node index.
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex3D(int const k) {
    return linearIndex3DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * ((k % 4) / 2), (num1dNodes - 1) * (k / 4));
  }

  /**
   * @brief Face-local index of the node (qa, qb).
   * @return qa + qb*num1dNodes.
   * @see docs/design.md, "Hexahedron local numbering".
   */
  PROXY_HOST_DEVICE
  constexpr static int linearIndex2DVal(const int qa, const int qb) { return qa + qb * num1dNodes; }

  /**
   * @brief Face-local index of the node at a vertex of a quadrilateral face.
   * @param[in] k Vertex index in [0, 3]; bits 0 and 1 give the side (0 =
   * minus, 1 = plus) along the first and second face axes.
   * @return The face-local node index.
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex2D(int const k) {
    return linearIndex2DVal((num1dNodes - 1) * (k % 2), (num1dNodes - 1) * (k / 2));
  }

  /**
   * @brief Weight of an endpoint of [-1, 1] in the linear interpolation at the
   * 1D node q.
   * @param[in] q 1D node index.
   * @param[in] k Endpoint: 0 for -1, 1 for +1.
   * @return (1 - xi_q)/2 for k = 0, (1 + xi_q)/2 for k = 1, where xi_q is the
   * parent coordinate of node q.
   */
  constexpr static real_t interpolationCoord(const int q, const int k) {
    const real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0);
    return k == 0 ? (1.0 - alpha) : alpha;
  }

  /**
   * @brief Derivative of the 1D basis function q at the 1D node p, with respect
   * to the parent coordinate in [-1, 1].
   *
   * Valid for every p in [0, num1dNodes): the half p > halfNodes is obtained by
   * symmetry from GL_BASIS::gradientAt().
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
   * @brief One 1D factor of the derivative of a trilinear vertex shape
   * function at a quadrature point.
   *
   * The derivative along parent axis @p dir of the shape function of the vertex
   * with sides (k0, k1, k2) is the product over i = 0, 1, 2 of
   * jacobianCoefficient1D(q_i, i, k_i, dir).
   * @param[in] q 1D quadrature point index along axis @p i.
   * @param[in] i Parent axis of this factor, in [0, 2].
   * @param[in] k Vertex side along axis @p i: 0 for -1, 1 for +1.
   * @param[in] dir Parent axis of the derivative, in [0, 2].
   * @return -1/2 or 1/2 when i == dir, interpolationCoord(q, k) otherwise.
   */
  PROXY_HOST_DEVICE
  constexpr static real_t jacobianCoefficient1D(const int q, const int i, const int k, const int dir) {
    if (i == dir) {
      return k == 0 ? -1.0 / 2.0 : 1.0 / 2.0;
    } else {
      return interpolationCoord(q, k);
    }
  }

  /// 1D Gauss-Lobatto quadrature weight of node q, for the interval [-1, 1].
  PROXY_HOST_DEVICE
  constexpr static real_t quadratureWeight(const int q) { return GL_BASIS::weight(q); }

  /// Number of quadrature points per element.
  PROXY_HOST_DEVICE static constexpr int getNumQuadraturePoints() { return numQuadraturePoints; }
  /// Number of nodes per element.
  PROXY_HOST_DEVICE static constexpr int getNumSupportPoints() { return numNodes; }
  /// Number of support points per element.
  PROXY_HOST_DEVICE static constexpr int getMaxSupportPoints() { return maxSupportPoints; }

  /**
   * @brief Adds the Jacobian of the bilinear map from [-1, 1]^2 to a
   * quadrilateral face, at the face quadrature point (qa, qb).
   * @param[in] X Coordinates of the 4 face vertices, X[vertex][axis], with the
   * vertex order of meshIndexToLinearIndex2D().
   * @param[in,out] J J[i][j] += d x_i / d xi_j; the caller zeroes it first.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation2d(int const qa, int const qb, real_t const (&X)[4][3], real_t (&J)[3][2]);

  /**
   * @brief Adds the Jacobian of the trilinear map from [-1, 1]^3 to the
   * hexahedron, at the quadrature point (qa, qb, qc).
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis], with the
   * vertex order of meshIndexToLinearIndex3D().
   * @param[in,out] J J[i][j] += d x_i / d xi_j; the caller zeroes it first.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                     real_t (&J)[3][3]);

  /**
   * @brief Computes the diagonal of the element mass matrix for a unit
   * coefficient: |det J| times the quadrature weight, at each node.
   * @param[in] X Coordinates of the 8 vertices, as in jacobianTransformation().
   * @param[in] func Called as func(q, value) once per node, q being the
   * element-local node index.
   */
  template <typename FUNC>
  PROXY_HOST_DEVICE static void computeMassTerm(float const (&X)[8][3], FUNC &&func);

  /**
   * @brief Surface measure times quadrature weight at a face node:
   * sqrt(det(J^T J)) * w, J being the 3x2 face Jacobian.
   * @param[in] q Face-local node index.
   * @param[in] X Coordinates of the 4 face vertices, as in
   * jacobianTransformation2d().
   */
  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]);

  /**
   * @brief Computes the Jacobian and the stiffness metric at the quadrature
   * point (qa, qb, qc).
   * @param[in] X Coordinates of the 8 vertices, as in jacobianTransformation().
   * @param[in,out] J Jacobian, accumulated as in jacobianTransformation(); the
   * caller zeroes it first.
   * @param[out] B det(J) * J^-1 J^-T in Voigt storage (docs/design.md,
   * "Symmetric 3x3 matrices"). det(J) is signed.
   */
  PROXY_HOST_DEVICE
  static void computeBMatrix(int const qa, int const qb, int const qc, real_t const (&X)[8][3], real_t (&J)[3][3],
                             real_t (&B)[6]);

  /**
   * @brief Enumerates the entries of the element stiffness matrix
   * K(i, j) = sum over quadrature points of w * grad(phi_i) . B grad(phi_j),
   * where grad is the gradient in parent coordinates, w the 3D quadrature
   * weight and B comes from computeBMatrix().
   * @param[in] X Coordinates of the 8 vertices, as in jacobianTransformation().
   * @param[in] func1 Called as func1(qa, qb, qc) at each quadrature point,
   * before the func2 calls of that point.
   * @param[in] func2 Called as func2(i, j, value), i and j being element-local
   * node indices; the calls for the same (i, j) must be summed.
   */
  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffnessTerm(float const (&X)[8][3], FUNC1 &&func1, FUNC2 &&func2);

  /**
   * @brief Adds K(alpha) p_local to f_local by sum factorization, without
   * forming the element stiffness matrix.
   *
   * K(alpha) is the stiffness matrix of computeStiffnessTerm() with each
   * quadrature point weighted by alpha.
   * @param[in] X Coordinates of the 8 vertices, as in jacobianTransformation().
   * @param[in] p_local Nodal values, indexed by element-local node index.
   * @param[in,out] f_local Result, accumulated (not zeroed).
   * @param[in] get_alpha Called as get_alpha(qa, qb, qc); returns the
   * coefficient at that quadrature point.
   */
  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(float const (&X)[8][3], real_t const (&p_local)[numNodes],
                                                            real_t (&f_local)[numNodes], FUNC_ALPHA &&get_alpha);

  /**
   * @brief Stiffness contributions of the quadrature point (qa, qb, qc); one
   * step of computeStiffnessTerm().
   * @param[in] B Metric at this point, from computeBMatrix().
   * @param[in] func1 See computeStiffnessTerm().
   * @param[in] func2 See computeStiffnessTerm().
   */
  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiBGradPhi(real_t const (&B)[6], FUNC1 &&func1, FUNC2 &&func2);

 protected:
  /// Protected and non-virtual: the class is never used through a base
  /// pointer, and a vtable would be unusable on device.
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
  // B = J^T J, 2x2 Voigt storage (B00, B11, B01).
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

  // B = J^T J / det(J), then inverted in place into det(J) J^-1 J^-T.
  B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) * invDetJ;
  B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) * invDetJ;
  B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]) * invDetJ;
  B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]) * invDetJ;
  B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]) * invDetJ;
  B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) * invDetJ;

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
      const real_t w0 = w * gia * gja;
      func2(ibc, jbc, w0 * B[0]);
      const real_t w1 = w * gib * gjb;
      func2(aic, ajc, w1 * B[1]);
      const real_t w2 = w * gic * gjc;
      func2(abi, abj, w2 * B[2]);
      // B is symmetric: each off-diagonal term contributes to (i, j) and (j, i).
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
  // Weighted fluxes at each quadrature point q:
  // (G_xi, G_eta, G_zeta)[q] = w_q * alpha_q * B_q * (parent gradient of u at q).
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

  // Apply the transposed 1D derivative matrix along each axis:
  // v[node] += sum over the quadrature points on the node's three lines.
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
