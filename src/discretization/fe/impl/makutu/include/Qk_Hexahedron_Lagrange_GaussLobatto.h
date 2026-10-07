/**
 * @file Qk_Hexahedron_Lagrange_GaussLobatto.h
 * @brief Element kernels of the makutu back-end: Qk hexahedra with
 * Gauss-Lobatto-Legendre nodes, and the compile-time loop helpers
 * for_constexpr() and triple_loop().
 */

#ifndef _QkHEXAHEDRON_HPP_
#define _QkHEXAHEDRON_HPP_

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
 * @brief Spectral-element kernels of a Qk hexahedron whose nodes, and
 * quadrature points, are the tensor-product Gauss-Lobatto-Legendre points.
 *
 * A stateless collection of static host/device functions called from the
 * solver kernels, one element at a time: local index maps, Jacobian of the map
 * from the parent cube [-1, 1]^3, shape function gradients, and the element
 * mass, damping, stiffness and interface-flux terms. The mass matrix is
 * diagonal because nodes and quadrature points coincide. Select a degree with
 * Qk_Hexahedron_Lagrange_GaussLobatto_Selector.
 * @tparam GL_BASIS 1D Lagrange basis, with the interface described in
 * docs/design.md, "1D Lagrange bases".
 * @see docs/design.md, "Hexahedron local numbering" and "Element geometry".
 */
template <typename GL_BASIS>
class Qk_Hexahedron_Lagrange_GaussLobatto {
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

  /// A 3x3 matrix in single precision, used for the Jacobian and its inverse.
  struct JacobianType {
    float data[3][3];  ///< Matrix entries, data[row][column].
  };

  /// Tag: this back-end provides the z-deformed sum-factorization kernels.
  struct ZDeformedSumFact {};

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
   * @see docs/design.md, "Element geometry".
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

  /// Non-virtual although the class declares virtual functions (see
  /// docs/design-red-flags.md).
  PROXY_HOST_DEVICE
  ~Qk_Hexahedron_Lagrange_GaussLobatto() = default;

  /// Number of quadrature points per element.
  PROXY_HOST_DEVICE
  virtual int getNumQuadraturePoints() { return numQuadraturePoints; }

  /// Number of nodes per element.
  PROXY_HOST_DEVICE
  virtual int getNumSupportPoints() { return numNodes; }

  /// Number of support points per element.
  PROXY_HOST_DEVICE
  virtual int getMaxSupportPoints() const { return maxSupportPoints; }

  /**
   * @brief Values of the numNodes shape functions at a point of the parent
   * cube.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[out] N N[l] is the shape function of element-local node l.
   */
  PROXY_HOST_DEVICE
  static void calcN(double const (&coords)[3], double (&N)[numNodes]) { GL_BASIS::TensorProduct3D::value(coords, N); }

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

  /**
   * @brief Values of the numNodes shape functions at the quadrature point q:
   * 1 for node q, 0 elsewhere, since nodes and quadrature points coincide.
   * @param[in] q Element-local quadrature point index.
   * @param[out] N N[l] is the shape function of element-local node l.
   */

  PROXY_HOST_DEVICE
  static void calcN(int const q, real_t (&N)[numNodes]) {
    for (int a = 0; a < numNodes; ++a) {
      N[a] = 0;
    }
    N[q] = 1.0;
  }

  /**
   * @brief Physical gradients of the numNodes shape functions at the
   * quadrature point q.
   *
   * The element geometry is the trilinear map of the 8 vertex nodes taken from
   * @p X; the other nodes of @p X are ignored.
   * @param[in] q Element-local quadrature point index.
   * @param[in] X Coordinates of every element node, X[node][axis], indexed by
   * element-local node index.
   * @param[out] gradN gradN[l][i] = d phi_l / d x_i.
   * @return det(J), signed.
   */

  PROXY_HOST_DEVICE
  static real_t calcGradN(int const q, real_t const (&X)[numNodes][3], real_t (&gradN)[numNodes][3]);
  /**
   * @brief Physical gradients of the numNodes shape functions at a point of the
   * parent cube, for the isoparametric map defined by all element nodes.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[in] X Coordinates of every element node, X[node][axis], indexed by
   * element-local node index.
   * @param[out] gradN gradN[l][i] = d phi_l / d x_i.
   * @return det(J), signed.
   */

  PROXY_HOST_DEVICE
  static real_t calcGradN(real_t const (&coords)[3], real_t const (&X)[numNodes][3], real_t (&gradN)[numNodes][3]);

  /**
   * @brief Physical gradients of the numNodes shape functions at the
   * quadrature point q, for the trilinear map of the 8 vertices.
   * @param[in] q Element-local quadrature point index.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[out] gradN gradN[l][i] = d phi_l / d x_i.
   * @return det(J), signed.
   * @see docs/design.md, "Element geometry".
   */

  PROXY_HOST_DEVICE
  static real_t calcGradNWithCorners(int const q, real_t const (&X)[8][3], real_t (&gradN)[numNodes][3]);
  /**
   * @brief Physical gradients of the numNodes shape functions at a point of the
   * parent cube, for the trilinear map of the 8 vertices.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[out] gradN gradN[l][i] = d phi_l / d x_i.
   * @return det(J), signed.
   * @see docs/design.md, "Element geometry".
   */

  PROXY_HOST_DEVICE
  static real_t calcGradNWithCorners(real_t const (&coords)[3], real_t const (&X)[8][3], real_t (&gradN)[numNodes][3]);

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
   * @brief Inverse Jacobian of the trilinear map at the quadrature point
   * (qa, qb, qc).
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in,out] J Must be zero on entry; holds J^-1 on return,
   * J^-1[r][i] = d xi_r / d x_i.
   * @return det(J), signed.
   * @see docs/design.md, "Element geometry".
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                          real_t (&J)[3][3]) {
    jacobianTransformation(qa, qb, qc, X, J);
    return invert3x3(J);
  }

  /**
   * @brief Same as the (qa, qb, qc) overload, for the element-local quadrature
   * point index q.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in,out] J Must be zero on entry; holds J^-1 on return.
   * @return det(J), signed.
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianTransformation(int const q, real_t const (&X)[8][3], real_t (&J)[3][3]) {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
    return invJacobianTransformation(qa, qb, qc, X, J);
  }

  /**
   * @brief Adds the symmetric gradient of a nodal vector field at the
   * quadrature point q.
   * @param[in] q Element-local quadrature point index.
   * @param[in] invJ Inverse Jacobian at q, invJ[r][i] = d xi_r / d x_i.
   * @param[in] var Nodal field, var[node][component].
   * @param[in,out] grad Accumulated (not zeroed). In the Voigt order of
   * docs/design.md, "Symmetric 3x3 matrices", entry (a, b) receives
   * d var_a / d x_b + d var_b / d x_a for a != b (engineering shear strain)
   * and d var_a / d x_a for a == b.
   */
  PROXY_HOST_DEVICE
  static void symmetricGradient(int const q, real_t const (&invJ)[3][3], real_t const (&var)[numNodes][3],
                                real_t (&grad)[6]);

  /**
   * @brief Adds the gradient of a nodal vector field at the quadrature point q:
   * grad[i][j] += sum over nodes a of d N_a / d x_j * var[a][i].
   * @param[in] q Element-local quadrature point index.
   * @param[in] invJ Inverse Jacobian at q, invJ[r][i] = d xi_r / d x_i.
   * @param[in] var Nodal field, var[node][component].
   * @param[in,out] grad Accumulated (not zeroed).
   */
  PROXY_HOST_DEVICE
  static void gradient(int const q, real_t const (&invJ)[3][3], real_t const (&var)[numNodes][3], real_t (&grad)[3][3]);

  /**
   * @brief Adds the Jacobian of the trilinear map from [-1, 1]^3 to the
   * hexahedron, at the quadrature point (qa, qb, qc).
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in,out] J J[i][j] += d x_i / d xi_j; the caller zeroes it first.
   * @see docs/design.md, "Element geometry".
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, int const qc, real_t const (&X)[8][3],
                                     real_t (&J)[3][3]);

  /**
   * @brief Adds the Jacobian of the isoparametric map defined by all element
   * nodes, at a point of the parent cube.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[in] X Coordinates of every element node, X[node][axis], indexed by
   * element-local node index.
   * @param[in,out] J J[i][j] += d x_i / d xi_j; the caller zeroes it first.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation(real_t const (&coords)[3], real_t const (&X)[numNodes][3], real_t (&J)[3][3]);

  /**
   * @brief Adds the Jacobian of the trilinear map of the 8 vertices, at a point
   * of the parent cube.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in,out] J J[i][j] += d x_i / d xi_j; the caller zeroes it first.
   * @see docs/design.md, "Element geometry".
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformationWithCorners(real_t const (&coords)[3], real_t const (&X)[8][3], real_t (&J)[3][3]);

  /**
   * @brief Image of a point of the unit cube by the trilinear map of the 8
   * vertices.
   * @param[in] alpha Coordinate along the first parent axis, in [0, 1].
   * @param[in] beta Coordinate along the second parent axis, in [0, 1].
   * @param[in] gamma Coordinate along the third parent axis, in [0, 1].
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[out] coords Physical coordinates of the point.
   * @see docs/design.md, "Element geometry".
   */
  PROXY_HOST_DEVICE
  static void trilinearInterp(real_t const alpha, real_t const beta, real_t const gamma, real_t const (&X)[8][3],
                              real_t (&coords)[3]);

  /**
   * @brief Physical coordinates of every element node, by trilinear
   * interpolation of the 8 vertices.
   * @param[in] Xmesh Coordinates of the 8 vertices, Xmesh[vertex][axis].
   * @param[out] X X[node][axis], indexed by element-local node index.
   */
  PROXY_HOST_DEVICE
  static void computeLocalCoords(real_t const (&Xmesh)[8][3], real_t (&X)[numNodes][3]);

  /**
   * @brief Computes the diagonal of the element mass matrix for a unit
   * coefficient: |det J| times the quadrature weight, at each node.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in] func Called as func(q, value) once per node, q being the
   * element-local node index.
   */
  template <typename FUNC>
  PROXY_HOST_DEVICE static void computeMassTerm(float const (&X)[8][3], FUNC &&func);

  /**
   * @brief Surface measure times quadrature weight at a face node:
   * sqrt(det(J^T J)) * w, J being the 3x2 face Jacobian. This is the diagonal
   * entry of the face mass matrix.
   * @param[in] q Face-local node index.
   * @param[in] X Coordinates of the 4 face vertices, as in
   * jacobianTransformation2d().
   */
  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]);

  /**
   * @brief Enumerates the nonzero face quadrature terms of
   * phi_j (grad(phi_i) . n) at the face quadrature point j = (qa, qb).
   *
   * Each value is w * dS * grad(phi_i) . n at the point, with w the 2D
   * quadrature weight and dS the surface measure. Summing value * p_i over the
   * nodes of both callbacks therefore gives w * dS * grad(p) . n at the point.
   * The contraction with @p kNormal is done here, once per point, rather than
   * per physical direction in the callbacks.
   * @param[in] qa 1D index of the point along the first face axis.
   * @param[in] qb 1D index of the point along the second face axis.
   * @param[in] kDir Normal parent axis of the face, in [0, 2].
   * @param[in] kQFixed 1D index of the face along the normal axis: 0 or
   * num1dNodes - 1.
   * @param[in] kX Coordinates of the 4 face vertices, as in
   * jacobianTransformation2d().
   * @param[in] invJ3D Inverse volume Jacobian at the point,
   * invJ3D[r][i] = d xi_r / d x_i.
   * @param[in] kNormal Normal vector to contract with, outward for the element
   * of @p invJ3D. Not normalized by the function.
   * @param[in] func Called as func(i, j, value) for the num1dNodes nodes along
   * each face axis through the point, i and j being face-local node indices.
   * @param[in] funcNormal Called as funcNormal(m, j, value) for the nodes on the
   * line through the point along the normal axis, m in [0, num1dNodes) being
   * the 1D index along that axis (m = kQFixed is the face node itself).
   * @see docs/design.md, "Hexahedron local numbering".
   */
  template <typename FUNC, typename FUNC_NORMAL>
  PROXY_HOST_DEVICE static void computeGradPhiPhiAt(int const qa, int const qb, int const kDir, int const kQFixed,
                                                    real_t const (&kX)[4][3], real_t const (&invJ3D)[3][3],
                                                    real_t const (&kNormal)[3], FUNC &&func, FUNC_NORMAL &&funcNormal);

  /**
   * @brief computeGradPhiPhiAt() at the face quadrature point q of face
   * @p kFaceId, with the volume Jacobian computed from the element vertices.
   * @param[in] q Face-local node index, in [0, numNodesPerFace).
   * @param[in] kX Coordinates of the 4 face vertices, as in
   * jacobianTransformation2d().
   * @param[in] X8 Coordinates of the 8 element vertices, X8[vertex][axis].
   * @param[in] kFaceId Face of the element, a model::CubicFace value in [0, 5].
   * @param[in] kNormal See computeGradPhiPhiAt().
   * @param[in] func See computeGradPhiPhiAt(); always called with j == q.
   * @param[in] funcNormal See computeGradPhiPhiAt(); always called with
   * j == q.
   * @warning The volume Jacobian is evaluated at the parent point
   * (qa, qb, kQFixed) whatever the face axis (see docs/design-red-flags.md).
   * @see docs/design.md, "Hexahedron local numbering".
   */
  template <typename FUNC, typename FUNC_NORMAL>
  PROXY_HOST_DEVICE static void computeInterfaceFluxTermAt(int const q, real_t const (&kX)[4][3],
                                                           real_t const (&X8)[8][3], int const kFaceId,
                                                           real_t const (&kNormal)[3], FUNC &&func,
                                                           FUNC_NORMAL &&funcNormal);

  /**
   * @brief Computes the Jacobian and the stiffness metric at the quadrature
   * point (qa, qb, qc).
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
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
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
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
   * quadrature point weighted by alpha. Cost O(num1dNodes^4) per element
   * instead of O(num1dNodes^5) for computeStiffnessTerm().
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in] p_local Nodal values, indexed by element-local node index.
   * @param[in,out] f_local Result, accumulated (not zeroed).
   * @param[in] get_alpha Called as get_alpha(qa, qb, qc); returns the
   * coefficient at that quadrature point.
   */
  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFact(float const (&X)[8][3], real_t const (&p_local)[numNodes],
                                                            real_t (&f_local)[numNodes], FUNC_ALPHA &&get_alpha);

  /**
   * @brief Inverse Jacobian at (qa, qb, qc) of an element deformed along z only.
   *
   * The element edges along xi and eta stay parallel to x and y, so J is lower
   * triangular with constant J[0][0] and J[1][1]; only its z row varies.
   * @param[in] J00 dx/dxi, half the element size along x.
   * @param[in] J11 dy/deta, half the element size along y.
   * @param[in] Z z coordinate of the 8 vertices.
   * @param[out] invJ invJ[r][i] = d xi_r / d x_i.
   * @return det(J).
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianZDeformed(int const qa, int const qb, int const qc, real_t const J00, real_t const J11,
                                     real_t const (&Z)[8], real_t (&invJ)[3][3]);

  /**
   * @brief Metric B = det(J) J^-1 J^-T of an element deformed along z only.
   *
   * The element edges along xi and eta stay parallel to x and y, so x depends
   * on xi only and y on eta only: J is lower triangular with constant
   * J[0][0] and J[1][1], and only the z row varies over the element. B[5]
   * (xy) is then zero. Same result as computeBMatrix() on such an element.
   * @param[in] J00 dx/dxi, half the element size along x.
   * @param[in] J11 dy/deta, half the element size along y.
   * @param[in] Z z coordinate of the 8 vertices.
   * @param[out] B Voigt order [xx, yy, zz, yz, xz, xy].
   */
  PROXY_HOST_DEVICE
  static void computeBMatrixZDeformed(int const qa, int const qb, int const qc, real_t const J00, real_t const J11,
                                      real_t const (&Z)[8], real_t (&B)[6]);

  /**
   * @brief computeStiffnessTermSumFact() for an element deformed along z only, without the
   * flux arrays.
   *
   * The metric comes from computeBMatrixZDeformed() instead of the full Jacobian. Each quadrature point adds its three fluxes to @p v_local as soon as they are
   * computed, so only @p u_local and @p v_local stay live: this keeps a
   * one-thread-per-element kernel within its registers. @p v_local is
   * overwritten and must not alias @p u_local.
   * @param[in] J00 dx/dxi, half the element size along x.
   * @param[in] J11 dy/deta, half the element size along y.
   * @param[in] Z z coordinate of the 8 vertices.
   * @param[in] u_local Nodal values, indexed by element-local node index.
   * @param[out] v_local Result; must not alias @p u_local.
   * @param[in] get_alpha Called as get_alpha(q); returns the coefficient at the
   * quadrature point q.
   */
  template <typename FUNC_ALPHA>
  PROXY_HOST_DEVICE static void computeStiffnessTermSumFactZDeformedScatter(real_t const J00, real_t const J11,
                                                                            real_t const (&Z)[8],
                                                                            real_t const (&u_local)[numNodes],
                                                                            real_t (&v_local)[numNodes],
                                                                            FUNC_ALPHA &&get_alpha);

  /**
   * @brief Stiffness contributions of the quadrature point (qa, qb, qc); one
   * step of computeStiffnessTerm().
   * @param[in] B Metric at this point, from computeBMatrix().
   * @param[in] func1 See computeStiffnessTerm().
   * @param[in] func2 See computeStiffnessTerm().
   */
  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiBGradPhi(real_t const (&B)[6], FUNC1 &&func1, FUNC2 &&func2);

  /**
   * @brief Stiffness contributions of the quadrature point (qa, qb, qc), per
   * pair of parent axes, for a constitutive law applied by the caller.
   * @param[in,out] J On entry the Jacobian at the point; inverted in place,
   * J^-1[r][i] = d xi_r / d x_i.
   * @param[in] func1 Called once as func1(qa, qb, qc, invJ) before the func2
   * calls, invJ being J^-1.
   * @param[in] func2 Called as func2(i, j, value, p, r) with
   * value = w * det(J) * d phi_i / d xi_p * d phi_j / d xi_r, i and j being
   * element-local node indices, w the 3D quadrature weight and det(J) signed.
   * The calls for the same (i, j) must be summed.
   */
  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiGradPhi(JacobianType &J, FUNC1 &&func1, FUNC2 &&func2);

  /**
   * @brief Runs computeGradPhiGradPhi() at every quadrature point of the
   * element.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in] func1 See computeGradPhiGradPhi().
   * @param[in] func2 See computeGradPhiGradPhi().
   */
  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffNessTermwithJac(float const (&X)[8][3], FUNC1 &&func1, FUNC2 &&func2);

  /**
   * @brief Adds the elastic stiffness term of a nodal displacement to f_local
   * by sum factorization, in O(num1dNodes^4) per element; the constitutive law
   * is supplied by the caller.
   *
   * At each quadrature point q the kernel computes the parent-coordinate
   * gradient of the displacement and asks @p func1 for a flux; it then adds,
   * for every node l and component f, the sum over q and p of
   * w * det(J) * d phi_l / d xi_p * flux_q[p][f] to f_local[f][l], where w is
   * the 3D quadrature weight and det(J) is signed.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in] u_local Displacement, u_local[component][node], indexed by
   * element-local node index.
   * @param[in,out] f_local Result, f_local[component][node], accumulated (not
   * zeroed).
   * @param[in] func1 Called once per quadrature point as
   * func1(qa, qb, qc, J_inv, grad_u_ref, flux), with J_inv[r][i] =
   * d xi_r / d x_i, grad_u_ref[r][s] = d u_s / d xi_r, and flux[p][f], zeroed
   * on entry, to be filled.
   */
  template <typename FUNC1>
  PROXY_HOST_DEVICE static void computeElasticStiffnessSumFact(float const (&X)[8][3],
                                                               real_t const (&u_local)[3][numNodes],
                                                               real_t (&f_local)[3][numNodes], FUNC1 &&func1);

  /**
   * @brief computeElasticStiffnessSumFact() computed by a Kokkos team on one
   * element, with caller-provided (team scratch) buffers.
   *
   * Must be called by every thread of the team. Unlike the single-thread
   * overload, @p f_local is overwritten, not accumulated, and may alias
   * @p u_local: a team barrier separates the last read of @p u_local from the
   * first write of @p f_local. The function does not end with a barrier.
   * @tparam TEAM_MEMBER Kokkos team member type.
   * @param[in] team Team handle.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[in] u_local Displacement, size 3*numNodes, at u_local[component *
   * numNodes + node].
   * @param[out] f_local Result, same layout as @p u_local.
   * @param[out] F Scratch of size 9*numNodes, at F[(p * 3 + f) * numNodes + q].
   * @param[in] func1 Same as in computeElasticStiffnessSumFact().
   * @param[in] basis_tab Optional table filled by basisTableEntry(), e.g. in
   * team scratch; the 1D basis values are then read from it instead of being
   * evaluated per thread from runtime indices. nullptr evaluates them.
   */
  template <typename TEAM_MEMBER, typename FUNC1>
  PROXY_HOST_DEVICE static void computeElasticStiffnessSumFactTeam(TEAM_MEMBER const &team, float const (&X)[8][3],
                                                                   real_t const *u_local, real_t *f_local, real_t *F,
                                                                   FUNC1 &&func1, real_t const *basis_tab = nullptr);

  /**
   * @brief Same as the vertex overload, for an element whose Jacobian is
   * constant (affine map), given precomputed.
   *
   * The other parameters and the synchronization rules are those of the vertex
   * overload.
   * @param[in] geom Size 10: J^-1 in row-major order (geom[r * 3 + i] =
   * d xi_r / d x_i), then det(J).
   */
  template <typename TEAM_MEMBER, typename FUNC1>
  PROXY_HOST_DEVICE static void computeElasticStiffnessSumFactTeam(TEAM_MEMBER const &team, real_t const *geom,
                                                                   real_t const *u_local, real_t *f_local, real_t *F,
                                                                   FUNC1 &&func1, real_t const *basis_tab = nullptr);

  /**
   * @brief Same as the vertex overload, for an element deformed along z only
   * (see invJacobianZDeformedBilinear()).
   *
   * The other parameters and the synchronization rules are those of the vertex
   * overload; @p basis_tab is required.
   * @param[in] J00 dx/dxi, half the element size along x.
   * @param[in] J11 dy/deta, half the element size along y.
   * @param[in] Z z coordinate of the 8 vertices.
   */
  template <typename TEAM_MEMBER, typename FUNC1>
  PROXY_HOST_DEVICE static void computeElasticStiffnessSumFactTeamZDeformed(
      TEAM_MEMBER const &team, real_t const J00, real_t const J11, real_t const (&Z)[8], real_t const *u_local,
      real_t *f_local, real_t *F, FUNC1 &&func1, real_t const *basis_tab);

  /// Size of the table read by the team kernels: D[q][p], then w[q], then alpha[q].
  constexpr static int kBasisTableSize = num1dNodes * num1dNodes + 2 * num1dNodes;

  /**
   * @brief Entry @p idx of the 1D basis table read by the team kernels.
   *
   * [0, n^2): basisGradientAt(q, p) at q * n + p; [n^2, n^2 + n): weight(q);
   * [n^2 + n, n^2 + 2n): interpolationCoord(q, 1), the position of node q in
   * [0, 1]. With n = num1dNodes.
   */
  PROXY_HOST_DEVICE
  static real_t basisTableEntry(int const idx);

  /**
   * @brief Body shared by the team overloads; @p geom_at gives the inverse
   * Jacobian and det(J) at each quadrature point.
   * @tparam kTab True to read the 1D basis from @p basis_tab.
   * @tparam GEOM Callable real_t(int qa, int qb, int qc, real_t const (&alpha)[3], real_t (&invJ)[3][3]),
   * returns det(J); alpha holds the interpolation coordinates of qa, qb, qc.
   */
  template <bool kTab, typename TEAM_MEMBER, typename GEOM, typename FUNC1>
  PROXY_HOST_DEVICE static void elasticStiffnessSumFactTeamImpl(TEAM_MEMBER const &team, GEOM &&geom_at,
                                                                real_t const *u_local, real_t *f_local, real_t *F,
                                                                FUNC1 &&func1, real_t const *basis_tab);

  /**
   * @brief Trilinear Jacobian at a point given by its interpolation coordinates.
   * @param[in] alpha Position of the point in [0, 1] along each parent axis.
   * @param[in] X Coordinates of the 8 vertices, X[vertex][axis].
   * @param[out] J J[i][j] = d x_i / d xi_j.
   */
  PROXY_HOST_DEVICE
  static void trilinearJacobian(real_t const (&alpha)[3], real_t const (&X)[8][3], real_t (&J)[3][3]);

  /**
   * @brief invJacobianZDeformed() at a point given by its interpolation
   * coordinates instead of its node indices.
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianZDeformed(real_t const (&alpha)[3], real_t const J00, real_t const J11,
                                     real_t const (&Z)[8], real_t (&invJ)[3][3]);

  /**
   * @brief Half differences of the vertex z along each parent axis, the input of
   * invJacobianZDeformedBilinear().
   * @param[in] Z z coordinate of the 8 vertices, vertex k = ka + 2 kb + 4 kc.
   * @param[out] dZ dZ[d][s0 + 2 s1]: half difference along axis d, (s0, s1) the sides on
   * the two other axes in increasing order.
   */
  PROXY_HOST_DEVICE
  static void zDeformedHalfDifferences(real_t const (&Z)[8], real_t (&dZ)[3][4]);

  /**
   * @brief invJacobianZDeformed() from zDeformedHalfDifferences(): each derivative of z is a
   * bilinear interpolation of the half differences, instead of a sum over the 8 vertices.
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianZDeformedBilinear(real_t const (&alpha)[3], real_t const J00, real_t const J11,
                                             real_t const (&dZ)[3][4], real_t (&invJ)[3][3]);

  /**
   * @brief Physical gradients of the numNodes shape functions at the
   * quadrature point q, from the inverse Jacobian at q.
   * @param[in] q Element-local quadrature point index.
   * @param[in] invJ Inverse Jacobian at q, invJ[r][i] = d xi_r / d x_i.
   * @param[out] gradN gradN[l][i] = d phi_l / d x_i; every entry is written.
   */
  PROXY_HOST_DEVICE
  static void applyTransformationToParentGradients(int const q, real_t const (&invJ)[3][3],
                                                   real_t (&gradN)[numNodes][3]);

  /**
   * @brief Physical gradients of the numNodes shape functions at a point of the
   * parent cube, from the inverse Jacobian at that point.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[in] invJ Inverse Jacobian at @p coords, invJ[r][i] = d xi_r / d x_i.
   * @param[out] gradN gradN[l][i] = d phi_l / d x_i; every entry is written.
   */
  PROXY_HOST_DEVICE
  static void applyTransformationToParentGradients(real_t const (&coords)[3], real_t const (&invJ)[3][3],
                                                   real_t (&gradN)[numNodes][3]);

 private:
  /// Distance between the first two 1D nodes of the parent interval. Unused.
  constexpr static real_t parentLength = GL_BASIS::parentSupportCoord(1) - GL_BASIS::parentSupportCoord(0);

  /// parentLength cubed. Unused.
  constexpr static real_t parentVolume = parentLength * parentLength * parentLength;
  /**
   * @brief Calls func(dNdXi, nodeIndex, params...) for every element node,
   * with the parent gradient of its shape function at a point of the parent
   * cube.
   * @param[in] coords Parent coordinates, each in [-1, 1].
   * @param[in] func Called with dNdXi[r] = d phi_nodeIndex / d xi_r.
   * @param[in,out] params Forwarded to @p func.
   */
  template <typename FUNC, typename... PARAMS>
  PROXY_HOST_DEVICE static void supportLoop(real_t const (&coords)[3], FUNC &&func, PARAMS &&...params);
  /**
   * @brief Calls func(dNdXi, nodeIndex, params...) for the 3*num1dNodes - 2
   * nodes whose parent gradient can be nonzero at the quadrature point
   * q = (qa, qb, qc).
   *
   * Because nodes and quadrature points coincide, the parent gradient of the
   * shape function of node (a, b, c) vanishes at q unless the node lies on one
   * of the three parent-axis lines through q. Every other node is skipped, so
   * a @p func that assigns rather than accumulates leaves their entries
   * untouched: the caller zeroes them.
   * @param[in] q Element-local quadrature point index.
   * @param[in] func Called with dNdXi[r] = d phi_nodeIndex / d xi_r at q.
   * @param[in,out] params Forwarded to @p func.
   */
  template <typename FUNC, typename... PARAMS>

  PROXY_HOST_DEVICE static void supportLoop(int const q, FUNC &&func, PARAMS &&...params);
};

/// @cond Doxygen_Suppress

template <typename GL_BASIS>
template <typename FUNC, typename... PARAMS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::supportLoop(real_t const (&coords)[3],
                                                                                  FUNC &&func, PARAMS &&...params) {
  for (int c = 0; c < num1dNodes; ++c) {
    for (int b = 0; b < num1dNodes; ++b) {
      for (int a = 0; a < num1dNodes; ++a) {
        real_t const dNdXi[3] = {static_cast<real_t>(GL_BASIS::gradient(a, coords[0]) * GL_BASIS::value(b, coords[1]) *
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

template <typename GL_BASIS>
template <typename FUNC, typename... PARAMS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::supportLoop(int const q, FUNC &&func,
                                                                                  PARAMS &&...params) {
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);

  // Line along the first axis. The node q itself is visited here, with its
  // three nonzero components; the other nodes of the line only have dNdXi[0].
  for (int a = 0; a < num1dNodes; ++a) {
    real_t const dNdXi[3] = {basisGradientAt(a, qa), (a == qa) ? basisGradientAt(qb, qb) : real_t(0),
                             (a == qa) ? basisGradientAt(qc, qc) : real_t(0)};
    int const nodeIndex = linearIndex3DVal(a, qb, qc);
    func(dNdXi, nodeIndex, std::forward<PARAMS>(params)...);
  }
  // Lines along the second and third axes, without the node q.
  for (int b = 0; b < num1dNodes; ++b) {
    if (b == qb) continue;
    real_t const dNdXi[3] = {real_t(0), basisGradientAt(b, qb), real_t(0)};
    int const nodeIndex = linearIndex3DVal(qa, b, qc);
    func(dNdXi, nodeIndex, std::forward<PARAMS>(params)...);
  }
  for (int c = 0; c < num1dNodes; ++c) {
    if (c == qc) continue;
    real_t const dNdXi[3] = {real_t(0), real_t(0), basisGradientAt(c, qc)};
    int const nodeIndex = linearIndex3DVal(qa, qb, c);
    func(dNdXi, nodeIndex, std::forward<PARAMS>(params)...);
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradN(int const q,
                                                                                  real_t const (&X)[numNodes][3],
                                                                                  real_t (&gradN)[numNodes][3]) {
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
  real_t Xmesh[8][3] = {{0}};
  for (int k = 0; k < 8; k++) {
    const int nodeIndex = meshIndexToLinearIndex3D(k);
    for (int i = 0; i < 3; i++) {
      Xmesh[k][i] = X[nodeIndex][i];
    }
  }
  real_t J[3][3] = {{0}};

  jacobianTransformation(qa, qb, qc, Xmesh, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(q, J, gradN);

  return detJ;
}
template <typename GL_BASIS>
PROXY_HOST_DEVICE

    real_t
    Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradN(real_t const (&coords)[3], real_t const (&X)[numNodes][3],
                                                             real_t (&gradN)[numNodes][3]) {
  real_t J[3][3] = {{0}};

  jacobianTransformation(coords, X, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(coords, J, gradN);

  return detJ;
}
template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradNWithCorners(
    int const q, real_t const (&X)[8][3], real_t (&gradN)[numNodes][3]) {
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);

  real_t J[3][3] = {{0}};

  jacobianTransformation(qa, qb, qc, X, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(q, J, gradN);

  return detJ;
}
template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradNWithCorners(
    real_t const (&coords)[3], real_t const (&X)[8][3], real_t (&gradN)[numNodes][3]) {
  real_t J[3][3] = {{0}};

  jacobianTransformationWithCorners(coords, X, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(coords, J, gradN);

  return detJ;
}

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

/// @endcond

/// Implementation of for_constexpr(); @p Is is 0 to N-1.
template <int N, typename F, int... Is>
constexpr void for_constexpr_impl(F &&f, std::integer_sequence<int, Is...>) {
  (f(std::integral_constant<int, Is>{}), ...);
}

/**
 * @brief Unrolled loop: calls f(std::integral_constant<int, I>{}) for I = 0 to
 * N-1, in increasing order, so that the index is a compile-time constant in
 * @p f.
 */
template <int N, typename F>
constexpr void for_constexpr(F &&f) {
  for_constexpr_impl<N>(std::forward<F>(f), std::make_integer_sequence<int, N>{});
}

/**
 * @brief Unrolled triple loop: calls lambda(I, J, K) with
 * std::integral_constant arguments for I in [0, BoundI), J in [0, BoundJ) and
 * K in [0, BoundK), K varying fastest.
 */

template <int BoundI, int BoundJ, int BoundK, typename Lambda>
constexpr void triple_loop(Lambda &&lambda) {
  for_constexpr<BoundI>(
      [&](auto I) { for_constexpr<BoundJ>([&](auto J) { for_constexpr<BoundK>([&](auto K) { lambda(I, J, K); }); }); });
}

/// @cond Doxygen_Suppress
template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformation(int const qa, int const qb,
                                                                                             int const qc,
                                                                                             real_t const (&X)[8][3],
                                                                                             real_t (&J)[3][3]) {
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
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformation(
    real_t const (&coords)[3], real_t const (&X)[numNodes][3], real_t (&J)[3][3]) {
  supportLoop(
      coords,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&X)[numNodes][3], real_t(&J)[3][3]) {
        real_t const *Xnode = X[nodeIndex];
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            J[i][j] = J[i][j] + dNdXi[j] * Xnode[i];
          }
        }
      },
      X, J);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformationWithCorners(
    real_t const (&coords)[3], real_t const (&X)[8][3], real_t (&J)[3][3]) {
  supportLoop(
      coords,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&X)[8][3], real_t(&J)[3][3]) {
        int qa, qb, qc;
        GL_BASIS::TensorProduct3D::multiIndex(nodeIndex, qa, qb, qc);
        real_t Xnode[3];
        real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(qa) + 1.0) / 2.0);
        real_t beta = static_cast<real_t>((GL_BASIS::parentSupportCoord(qb) + 1.0) / 2.0);
        real_t gamma = static_cast<real_t>((GL_BASIS::parentSupportCoord(qc) + 1.0) / 2.0);
        trilinearInterp(alpha, beta, gamma, X, Xnode);
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            J[i][j] = J[i][j] + dNdXi[j] * Xnode[i];
          }
        }
      },
      X, J);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::trilinearInterp(
    real_t const alpha, real_t const beta, real_t const gamma, real_t const (&X)[8][3], real_t (&coords)[3]) {
  for (int i = 0; i < 3; i++) {
    coords[i] = X[0][i] * (1.0 - alpha) * (1.0 - beta) * (1.0 - gamma) +
                X[1][i] * alpha * (1.0 - beta) * (1.0 - gamma) + X[2][i] * (1.0 - alpha) * beta * (1.0 - gamma) +
                X[3][i] * alpha * beta * (1.0 - gamma) + X[4][i] * (1.0 - alpha) * (1.0 - beta) * gamma +
                X[5][i] * alpha * (1.0 - beta) * gamma + X[6][i] * (1.0 - alpha) * beta * gamma +
                X[7][i] * alpha * beta * gamma;
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeLocalCoords(real_t const (&Xmesh)[8][3],
                                                                                         real_t (&X)[numNodes][3]) {
  int qa, qb, qc;
  for (int q = 0; q < numNodes; q++) {
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
    real_t alpha = static_cast<real_t>((GL_BASIS::parentSupportCoord(qa) + 1.0) / 2.0);
    real_t beta = static_cast<real_t>((GL_BASIS::parentSupportCoord(qb) + 1.0) / 2.0);
    real_t gamma = static_cast<real_t>((GL_BASIS::parentSupportCoord(qc) + 1.0) / 2.0);
    trilinearInterp(alpha, beta, gamma, Xmesh, X[q]);
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformation2d(int const qa,
                                                                                               int const qb,
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
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeMassTerm(float const (&X)[8][3],
                                                                                      FUNC &&func) {
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
template <typename FUNC, typename FUNC_NORMAL>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeGradPhiPhiAt(
    int const qa, int const qb, int const kDir, int const kQFixed, real_t const (&kX)[4][3],
    real_t const (&invJ3D)[3][3], real_t const (&kNormal)[3], FUNC &&func, FUNC_NORMAL &&funcNormal) {
  int ifa, ifb;
  switch (kDir) {
    case 0:
      ifa = 1;
      ifb = 2;
      break;
    case 1:
      ifa = 0;
      ifb = 2;
      break;
    default:
      ifa = 0;
      ifb = 1;
      break;
  }
  // Narrowed before the product: weight() returns double in several bases, and with runtime indices
  // there is no constant folding left to absorb it.
  const real_t kW2D = static_cast<real_t>(GL_BASIS::weight(qa)) * static_cast<real_t>(GL_BASIS::weight(qb));
  real_t B[3];
  real_t J[3][2] = {{0}};
  jacobianTransformation2d(qa, qb, kX, J);
  // B = J^T J, 2x2 Voigt storage (B00, B11, B01).
  B[0] = J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0];
  B[1] = J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1];
  B[2] = J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1];
  const real_t kDetJ = std::sqrt(std::abs(symDeterminant(B)));
  const real_t kVal = kW2D * kDetJ;
  const int kAbj = GL_BASIS::TensorProduct2D::linearIndex(qa, qb);

  // (J^-1 n) along each parent axis: d phi / d xi_r times these gives grad(phi) . n.
  const real_t kSa = invJ3D[ifa][0] * kNormal[0] + invJ3D[ifa][1] * kNormal[1] + invJ3D[ifa][2] * kNormal[2];
  const real_t kSb = invJ3D[ifb][0] * kNormal[0] + invJ3D[ifb][1] * kNormal[1] + invJ3D[ifb][2] * kNormal[2];
  const real_t kSd = invJ3D[kDir][0] * kNormal[0] + invJ3D[kDir][1] * kNormal[1] + invJ3D[kDir][2] * kNormal[2];

  for (int i = 0; i < num1dNodes; i++) {
    const int kIb = GL_BASIS::TensorProduct2D::linearIndex(i, qb);
    const int kAi = GL_BASIS::TensorProduct2D::linearIndex(qa, i);
    func(kIb, kAbj, kVal * kSa * basisGradientAt(i, qa));
    func(kAi, kAbj, kVal * kSb * basisGradientAt(i, qb));
    funcNormal(i, kAbj, kVal * kSd * basisGradientAt(i, kQFixed));
  }
}

template <typename GL_BASIS>
template <typename FUNC, typename FUNC_NORMAL>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeInterfaceFluxTermAt(
    int const q, real_t const (&kX)[4][3], real_t const (&X8)[8][3], int const kFaceId, real_t const (&kNormal)[3],
    FUNC &&func, FUNC_NORMAL &&funcNormal) {
  const int kDir = kFaceId / 2;
  const int kQFixed = (kFaceId % 2 == 0) ? 0 : num1dNodes - 1;
  int qa, qb;
  GL_BASIS::TensorProduct2D::multiIndex(q, qa, qb);
  real_t invJ3D[3][3] = {{0}};
  invJacobianTransformation(qa, qb, kQFixed, X8, invJ3D);
  computeGradPhiPhiAt(qa, qb, kDir, kQFixed, kX, invJ3D, kNormal, func, funcNormal);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeDampingTerm(int const q,
                                                                                           real_t const (&X)[4][3]) {
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
  return std::sqrt(std::abs(symDeterminant(B))) * w2D;
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeBMatrix(
    int const qa, int const qb, int const qc, real_t const (&X)[8][3], real_t (&J)[3][3], real_t (&B)[6]) {
  jacobianTransformation(qa, qb, qc, X, J);
  real_t const detJ = determinant(J);
  real_t const invDetJ = real_t(1) / detJ;

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
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeGradPhiBGradPhi(real_t const (&B)[6],
                                                                                             FUNC1 &&func1,
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
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeStiffnessTerm(float const (&X)[8][3],
                                                                                           FUNC1 &&func1,
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
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeStiffnessTermSumFact(
    float const (&X)[8][3], real_t const (&u_local)[numNodes], real_t (&v_local)[numNodes], FUNC_ALPHA &&get_alpha) {
  // Weighted parent-coordinate fluxes at each quadrature point q:
  // (G_xi, G_eta, G_zeta)[q] = w * alpha * B (parent gradient of u).
  real_t G_xi[numNodes] = {0};
  real_t G_eta[numNodes] = {0};
  real_t G_zeta[numNodes] = {0};

  // Passes 1 and 2: parent gradient of u, then the fluxes.
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

  // Pass 3: v[node] += sum over q and axis r of d phi_node / d xi_r (q) * G_r[q],
  // restricted to the quadrature points on the lines through the node.
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

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::basisTableEntry(int const idx) {
  constexpr int n = num1dNodes;
  if (idx < n * n) return basisGradientAt(idx / n, idx % n);
  if (idx < n * n + n) return static_cast<real_t>(GL_BASIS::weight(idx - n * n));
  return interpolationCoord(idx - n * n - n, 1);
}

/// Factor of the trilinear shape function derivative along direction j, for
/// parent axis d, vertex side kd (0 or 1) and interpolation coordinate a.
PROXY_HOST_DEVICE real_t qkTrilinearFactor(int const d, int const kd, int const j, real_t const a) {
  if (d == j) return kd == 0 ? real_t(-0.5) : real_t(0.5);
  return kd == 0 ? real_t(1) - a : a;
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::trilinearJacobian(real_t const (&alpha)[3],
                                                                                        real_t const (&X)[8][3],
                                                                                        real_t (&J)[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) J[i][j] = real_t(0);
  for (int k = 0; k < 8; k++) {
    const int ka = k % 2;
    const int kb = (k % 4) / 2;
    const int kc = k / 4;
    for (int j = 0; j < 3; j++) {
      real_t const c = qkTrilinearFactor(0, ka, j, alpha[0]) * qkTrilinearFactor(1, kb, j, alpha[1]) *
                       qkTrilinearFactor(2, kc, j, alpha[2]);
      for (int i = 0; i < 3; i++) J[i][j] += c * X[k][i];
    }
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::invJacobianZDeformed(
    int const qa, int const qb, int const qc, real_t const J00, real_t const J11, real_t const (&Z)[8],
    real_t (&invJ)[3][3]) {
  real_t const alpha[3] = {interpolationCoord(qa, 1), interpolationCoord(qb, 1), interpolationCoord(qc, 1)};
  return invJacobianZDeformed(alpha, J00, J11, Z, invJ);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::invJacobianZDeformed(
    real_t const (&alpha)[3], real_t const J00, real_t const J11, real_t const (&Z)[8], real_t (&invJ)[3][3]) {
  real_t J2[3] = {0, 0, 0};
  for (int k = 0; k < 8; k++) {
    const int ka = k % 2;
    const int kb = (k % 4) / 2;
    const int kc = k / 4;
    for (int j = 0; j < 3; j++) {
      J2[j] += qkTrilinearFactor(0, ka, j, alpha[0]) * qkTrilinearFactor(1, kb, j, alpha[1]) *
               qkTrilinearFactor(2, kc, j, alpha[2]) * Z[k];
    }
  }

  // J^-1 rows: (a, 0, 0), (0, b, 0), (c, d, e).
  real_t const a = real_t(1) / J00;
  real_t const b = real_t(1) / J11;
  real_t const e = real_t(1) / J2[2];
  invJ[0][0] = a;
  invJ[0][1] = real_t(0);
  invJ[0][2] = real_t(0);
  invJ[1][0] = real_t(0);
  invJ[1][1] = b;
  invJ[1][2] = real_t(0);
  invJ[2][0] = -a * J2[0] * e;
  invJ[2][1] = -b * J2[1] * e;
  invJ[2][2] = e;
  return J00 * J11 * J2[2];
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::zDeformedHalfDifferences(real_t const (&Z)[8],
                                                                                              real_t (&dZ)[3][4]) {
  for (int s1 = 0; s1 < 2; ++s1)
    for (int s0 = 0; s0 < 2; ++s0) {
      dZ[0][s0 + 2 * s1] = real_t(0.5) * (Z[1 + 2 * s0 + 4 * s1] - Z[2 * s0 + 4 * s1]);
      dZ[1][s0 + 2 * s1] = real_t(0.5) * (Z[s0 + 2 + 4 * s1] - Z[s0 + 4 * s1]);
      dZ[2][s0 + 2 * s1] = real_t(0.5) * (Z[s0 + 2 * s1 + 4] - Z[s0 + 2 * s1]);
    }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::invJacobianZDeformedBilinear(
    real_t const (&alpha)[3], real_t const J00, real_t const J11, real_t const (&dZ)[3][4], real_t (&invJ)[3][3]) {
  real_t const a0 = alpha[0], a1 = alpha[1], a2 = alpha[2];
  real_t const J20 = (real_t(1) - a2) * ((real_t(1) - a1) * dZ[0][0] + a1 * dZ[0][1]) +
                     a2 * ((real_t(1) - a1) * dZ[0][2] + a1 * dZ[0][3]);
  real_t const J21 = (real_t(1) - a2) * ((real_t(1) - a0) * dZ[1][0] + a0 * dZ[1][1]) +
                     a2 * ((real_t(1) - a0) * dZ[1][2] + a0 * dZ[1][3]);
  real_t const J22 = (real_t(1) - a1) * ((real_t(1) - a0) * dZ[2][0] + a0 * dZ[2][1]) +
                     a1 * ((real_t(1) - a0) * dZ[2][2] + a0 * dZ[2][3]);

  // J^-1 rows: (a, 0, 0), (0, b, 0), (c, d, e).
  real_t const a = real_t(1) / J00;
  real_t const b = real_t(1) / J11;
  real_t const e = real_t(1) / J22;
  invJ[0][0] = a;
  invJ[0][1] = real_t(0);
  invJ[0][2] = real_t(0);
  invJ[1][0] = real_t(0);
  invJ[1][1] = b;
  invJ[1][2] = real_t(0);
  invJ[2][0] = -a * J20 * e;
  invJ[2][1] = -b * J21 * e;
  invJ[2][2] = e;
  return J00 * J11 * J22;
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeBMatrixZDeformed(
    int const qa, int const qb, int const qc, real_t const J00, real_t const J11, real_t const (&Z)[8],
    real_t (&B)[6]) {
  real_t invJ[3][3];
  real_t const detJ = invJacobianZDeformed(qa, qb, qc, J00, J11, Z, invJ);
  real_t const a = invJ[0][0];
  real_t const b = invJ[1][1];
  real_t const c = invJ[2][0];
  real_t const d = invJ[2][1];
  real_t const e = invJ[2][2];

  B[0] = detJ * a * a;
  B[1] = detJ * b * b;
  B[2] = detJ * (c * c + d * d + e * e);
  B[3] = detJ * b * d;
  B[4] = detJ * a * c;
  B[5] = real_t(0);
}

template <typename GL_BASIS>
template <typename FUNC_ALPHA>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeStiffnessTermSumFactZDeformedScatter(
    real_t const J00, real_t const J11, real_t const (&Z)[8], real_t const (&u_local)[numNodes],
    real_t (&v_local)[numNodes], FUNC_ALPHA &&get_alpha) {
  for_constexpr<numNodes>([&](auto icn) { v_local[decltype(icn)::value] = 0; });

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

    real_t B[6];
    computeBMatrixZDeformed(qa, qb, qc, J00, J11, Z, B);

    real_t const scale = w * get_alpha(q);

    // B[5] (xy) is zero.
    real_t const g_xi = scale * (B[0] * dxi_q + B[4] * dzeta_q);
    real_t const g_eta = scale * (B[1] * deta_q + B[3] * dzeta_q);
    real_t const g_zeta = scale * (B[4] * dxi_q + B[3] * deta_q + B[2] * dzeta_q);

    // Transpose of the gradient above: node (i, qb, qc) receives D(i, qa) g_xi, and so on.
    for_constexpr<num1dNodes>([&](auto ici) {
      constexpr int i = decltype(ici)::value;
      v_local[GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc)] += basisGradientAt(i, qa) * g_xi;
      v_local[GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc)] += basisGradientAt(i, qb) * g_eta;
      v_local[GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i)] += basisGradientAt(i, qc) * g_zeta;
    });
  });
}

template <typename GL_BASIS>
template <typename FUNC1, typename FUNC2>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeStiffNessTermwithJac(
    float const (&X)[8][3], FUNC1 &&func1, FUNC2 &&func2) {
  triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
    constexpr int qa = decltype(icqa)::value;
    constexpr int qb = decltype(icqb)::value;
    constexpr int qc = decltype(icqc)::value;
    JacobianType J = {{0}};
    jacobianTransformation(qa, qb, qc, X, J.data);
    computeGradPhiGradPhi<qa, qb, qc>(J, func1, func2);
  });
}

template <typename GL_BASIS>
template <typename FUNC1>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeElasticStiffnessSumFact(
    float const (&X)[8][3], real_t const (&u_local)[3][numNodes], real_t (&f_local)[3][numNodes], FUNC1 &&func1) {
  // Scaled fluxes, one array per parent axis p, indexed [component f][q].
  real_t F_xi[3][numNodes] = {{0}};
  real_t F_eta[3][numNodes] = {{0}};
  real_t F_zeta[3][numNodes] = {{0}};

  // Passes 1 and 2: parent gradient of u, then the fluxes from func1.
  triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icqa, auto const icqb, auto const icqc) {
    constexpr int qa = decltype(icqa)::value;
    constexpr int qb = decltype(icqb)::value;
    constexpr int qc = decltype(icqc)::value;
    constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
    constexpr real_t w = GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);

    real_t grad_u_ref[3][3] = {{0}};
    for_constexpr<num1dNodes>([&](auto ici) {
      constexpr int i = decltype(ici)::value;
      constexpr int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
      constexpr int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
      constexpr int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
      const real_t gxi = basisGradientAt(i, qa);
      const real_t geta = basisGradientAt(i, qb);
      const real_t gzeta = basisGradientAt(i, qc);
      for (int s = 0; s < 3; ++s) {
        grad_u_ref[0][s] += gxi * u_local[s][ibc];
        grad_u_ref[1][s] += geta * u_local[s][aic];
        grad_u_ref[2][s] += gzeta * u_local[s][abi];
      }
    });

    JacobianType J = {{0}};
    jacobianTransformation(qa, qb, qc, X, J.data);
    real_t const detJ = invert3x3(J.data);
    const real_t scale = w * detJ;

    real_t flux[3][3] = {{0}};
    func1(qa, qb, qc, J.data, grad_u_ref, flux);

    F_xi[0][q] = scale * flux[0][0];
    F_xi[1][q] = scale * flux[0][1];
    F_xi[2][q] = scale * flux[0][2];
    F_eta[0][q] = scale * flux[1][0];
    F_eta[1][q] = scale * flux[1][1];
    F_eta[2][q] = scale * flux[1][2];
    F_zeta[0][q] = scale * flux[2][0];
    F_zeta[1][q] = scale * flux[2][1];
    F_zeta[2][q] = scale * flux[2][2];
  });

  // Pass 3: f_local[f][node] += sum over q and axis p of
  // d phi_node / d xi_p (q) * F_p[f][q]. The f loop is innermost so that the
  // basis derivative is computed once for the three components.
  triple_loop<num1dNodes, num1dNodes, num1dNodes>([&](auto const icia, auto const icib, auto const icic) {
    constexpr int ia = decltype(icia)::value;
    constexpr int ib = decltype(icib)::value;
    constexpr int ic = decltype(icic)::value;
    constexpr int node = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, ic);

    real_t v[3] = {0};
    for_constexpr<num1dNodes>([&](auto icqa) {
      constexpr int qa = decltype(icqa)::value;
      constexpr int q_xi = GL_BASIS::TensorProduct3D::linearIndex(qa, ib, ic);
      const real_t g = basisGradientAt(ia, qa);
      for (int f = 0; f < 3; ++f) v[f] += g * F_xi[f][q_xi];
    });
    for_constexpr<num1dNodes>([&](auto icqb) {
      constexpr int qb = decltype(icqb)::value;
      constexpr int q_eta = GL_BASIS::TensorProduct3D::linearIndex(ia, qb, ic);
      const real_t g = basisGradientAt(ib, qb);
      for (int f = 0; f < 3; ++f) v[f] += g * F_eta[f][q_eta];
    });
    for_constexpr<num1dNodes>([&](auto icqc) {
      constexpr int qc = decltype(icqc)::value;
      constexpr int q_zeta = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, qc);
      const real_t g = basisGradientAt(ic, qc);
      for (int f = 0; f < 3; ++f) v[f] += g * F_zeta[f][q_zeta];
    });
    for (int f = 0; f < 3; ++f) f_local[f][node] += v[f];
  });
}

template <typename GL_BASIS>
template <bool kTab, typename TEAM_MEMBER, typename GEOM, typename FUNC1>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::elasticStiffnessSumFactTeamImpl(
    TEAM_MEMBER const &team, GEOM &&geom_at, real_t const *u_local, real_t *f_local, real_t *F, FUNC1 &&func1,
    real_t const *basis_tab) {
  constexpr int n = num1dNodes;
  // With runtime indices, the basis functions compile to branches returning doubles; the table avoids both.
  auto const grad = [&](int const q, int const p) -> real_t {
    if constexpr (kTab) {
      return basis_tab[q * n + p];
    } else {
      return basisGradientAt(q, p);
    }
  };
  auto const weight1d = [&](int const q) -> real_t {
    if constexpr (kTab) {
      return basis_tab[n * n + q];
    } else {
      return static_cast<real_t>(GL_BASIS::weight(q));
    }
  };
  auto const alpha1d = [&](int const q) -> real_t {
    if constexpr (kTab) {
      return basis_tab[n * n + n + q];
    } else {
      return interpolationCoord(q, 1);
    }
  };

  // Pass 1+2: one thread per quadrature point. Reference gradients, then the
  // constitutive callback, then scale and store into the flux scratch.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numNodes), [&](const int q) {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t grad_u_ref[3][3] = {{0}};
    for (int i = 0; i < n; ++i) {
      int const ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
      int const aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
      int const abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
      real_t const gxi = grad(i, qa);
      real_t const geta = grad(i, qb);
      real_t const gzeta = grad(i, qc);
      for (int s = 0; s < 3; ++s) {
        grad_u_ref[0][s] += gxi * u_local[s * numNodes + ibc];
        grad_u_ref[1][s] += geta * u_local[s * numNodes + aic];
        grad_u_ref[2][s] += gzeta * u_local[s * numNodes + abi];
      }
    }

    real_t const alpha[3] = {alpha1d(qa), alpha1d(qb), alpha1d(qc)};
    real_t J_inv[3][3];
    real_t const detJ = geom_at(qa, qb, qc, alpha, J_inv);
    real_t const scale = weight1d(qa) * weight1d(qb) * weight1d(qc) * detJ;

    real_t flux[3][3] = {{0}};
    func1(qa, qb, qc, J_inv, grad_u_ref, flux);

    for (int p = 0; p < 3; ++p)
      for (int f = 0; f < 3; ++f) F[(p * 3 + f) * numNodes + q] = scale * flux[p][f];
  });
  team.team_barrier();

  // Pass 3: one thread per node, contract D^T with the stored fluxes. Each node
  // is owned by exactly one thread, so the accumulation needs no atomics.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numNodes), [&](const int node) {
    int ia, ib, ic;
    GL_BASIS::TensorProduct3D::multiIndex(node, ia, ib, ic);

    real_t v[3] = {0};
    for (int qa = 0; qa < n; ++qa) {
      int const q_xi = GL_BASIS::TensorProduct3D::linearIndex(qa, ib, ic);
      real_t const g = grad(ia, qa);
      for (int f = 0; f < 3; ++f) v[f] += g * F[(0 * 3 + f) * numNodes + q_xi];
    }
    for (int qb = 0; qb < n; ++qb) {
      int const q_eta = GL_BASIS::TensorProduct3D::linearIndex(ia, qb, ic);
      real_t const g = grad(ib, qb);
      for (int f = 0; f < 3; ++f) v[f] += g * F[(1 * 3 + f) * numNodes + q_eta];
    }
    for (int qc = 0; qc < n; ++qc) {
      int const q_zeta = GL_BASIS::TensorProduct3D::linearIndex(ia, ib, qc);
      real_t const g = grad(ic, qc);
      for (int f = 0; f < 3; ++f) v[f] += g * F[(2 * 3 + f) * numNodes + q_zeta];
    }
    for (int f = 0; f < 3; ++f) f_local[f * numNodes + node] = v[f];
  });
}

template <typename GL_BASIS>
template <typename TEAM_MEMBER, typename FUNC1>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeElasticStiffnessSumFactTeam(
    TEAM_MEMBER const &team, float const (&X)[8][3], real_t const *u_local, real_t *f_local, real_t *F, FUNC1 &&func1,
    real_t const *basis_tab) {
  auto const geom_at = [&](int, int, int, real_t const(&alpha)[3], real_t(&J_inv)[3][3]) -> real_t {
    trilinearJacobian(alpha, X, J_inv);
    return invert3x3(J_inv);
  };
  if (basis_tab != nullptr)
    elasticStiffnessSumFactTeamImpl<true>(team, geom_at, u_local, f_local, F, func1, basis_tab);
  else
    elasticStiffnessSumFactTeamImpl<false>(team, geom_at, u_local, f_local, F, func1, basis_tab);
}

template <typename GL_BASIS>
template <typename TEAM_MEMBER, typename FUNC1>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeElasticStiffnessSumFactTeam(
    TEAM_MEMBER const &team, real_t const *geom, real_t const *u_local, real_t *f_local, real_t *F, FUNC1 &&func1,
    real_t const *basis_tab) {
  // Geometry is element-wide: read it instead of rebuilding it per point.
  auto const geom_at = [&](int, int, int, real_t const(&)[3], real_t(&J_inv)[3][3]) -> real_t {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) J_inv[a][b] = geom[a * 3 + b];
    return geom[9];
  };
  if (basis_tab != nullptr)
    elasticStiffnessSumFactTeamImpl<true>(team, geom_at, u_local, f_local, F, func1, basis_tab);
  else
    elasticStiffnessSumFactTeamImpl<false>(team, geom_at, u_local, f_local, F, func1, basis_tab);
}

template <typename GL_BASIS>
template <typename TEAM_MEMBER, typename FUNC1>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeElasticStiffnessSumFactTeamZDeformed(
    TEAM_MEMBER const &team, real_t const J00, real_t const J11, real_t const (&Z)[8], real_t const *u_local,
    real_t *f_local, real_t *F, FUNC1 &&func1, real_t const *basis_tab) {
  real_t dZ[3][4];
  zDeformedHalfDifferences(Z, dZ);
  auto const geom_at = [&](int, int, int, real_t const(&alpha)[3], real_t(&J_inv)[3][3]) -> real_t {
    return invJacobianZDeformedBilinear(alpha, J00, J11, dZ, J_inv);
  };
  elasticStiffnessSumFactTeamImpl<true>(team, geom_at, u_local, f_local, F, func1, basis_tab);
}

template <typename GL_BASIS>
template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeGradPhiGradPhi(JacobianType &J,
                                                                                            FUNC1 &&func1,
                                                                                            FUNC2 &&func2) {
  real_t const detJ = invert3x3(J.data);
  const real_t w = static_cast<real_t>(GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc));
  func1(qa, qb, qc, J.data);
#pragma unroll 1
  for (int i = 0; i < num1dNodes; i++) {
    const int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
    const int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
    const int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
    const real_t gia = basisGradientAt(i, qa);
    const real_t gib = basisGradientAt(i, qb);
    const real_t gic = basisGradientAt(i, qc);
#pragma unroll 1
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
      // Each off-diagonal pair of parent axes contributes to (i, j) and (j, i).
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

template <typename GL_BASIS>

PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::applyTransformationToParentGradients(
    int const q, real_t const (&invJ)[3][3], real_t (&gradN)[numNodes][3]) {
  // The sparse supportLoop only visits the 3N-2 non-zero nodes; off-line
  // entries must be zeroed explicitly since this function assigns (not
  // accumulates) gradN[nodeIndex].
  for (int node = 0; node < numNodes; ++node) {
    gradN[node][0] = gradN[node][1] = gradN[node][2] = real_t(0);
  }
  supportLoop(
      q,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&invJ)[3][3], real_t(&gradN)[numNodes][3]) {
        // Unrolled by hand to reduce register pressure.
        gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
        gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
        gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];
      },
      invJ, gradN);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::applyTransformationToParentGradients(
    real_t const (&coords)[3], real_t const (&invJ)[3][3], real_t (&gradN)[numNodes][3]) {
  supportLoop(
      coords,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&invJ)[3][3], real_t(&gradN)[numNodes][3]) {
        gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
        gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
        gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];
      },
      invJ, gradN);
}

template <typename GL_BASIS>

PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::symmetricGradient(
    int const q, real_t const (&invJ)[3][3], real_t const (&var)[numNodes][3], real_t (&grad)[6]) {
  supportLoop(
      q,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&invJ)[3][3], real_t const(&var)[numNodes][3],
         real_t(&grad)[6]) {
        real_t gradN[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            gradN[i] = gradN[i] + dNdXi[j] * invJ[j][i];
          }
        }

        grad[0] = grad[0] + gradN[0] * var[nodeIndex][0];
        grad[1] = grad[1] + gradN[1] * var[nodeIndex][1];
        grad[2] = grad[2] + gradN[2] * var[nodeIndex][2];
        grad[3] = grad[3] + gradN[2] * var[nodeIndex][1] + gradN[1] * var[nodeIndex][2];
        grad[4] = grad[4] + gradN[2] * var[nodeIndex][0] + gradN[0] * var[nodeIndex][2];
        grad[5] = grad[5] + gradN[1] * var[nodeIndex][0] + gradN[0] * var[nodeIndex][1];
      },
      invJ, var, grad);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::gradient(int const q, real_t const (&invJ)[3][3],
                                                                               real_t const (&var)[numNodes][3],
                                                                               real_t (&grad)[3][3]) {
  supportLoop(
      q,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&invJ)[3][3], real_t const(&var)[numNodes][3],
         real_t(&grad)[3][3]) {
        for (int i = 0; i < 3; ++i) {
          real_t gradN = 0.0;
          ;
          for (int j = 0; j < 3; ++j) {
            gradN = gradN + dNdXi[j] * invJ[j][i];
          }
          for (int k = 0; k < 3; ++k) {
            grad[k][i] = grad[k][i] + gradN * var[nodeIndex][k];
          }
        }
      },
      invJ, var, grad);
}

/// @endcond

/// Q1 hexahedron (8 nodes).
using Q1_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis1>;

/// Q2 hexahedron (3^3 nodes).
using Q2_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis2>;

/// Q3 hexahedron (4^3 nodes).
using Q3_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis3GL>;

/// Q4 hexahedron (5^3 nodes).
using Q4_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis4GL>;

/// Q5 hexahedron (6^3 nodes).
using Q5_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis5GL>;

/// Q6 hexahedron (7^3 nodes).
using Q6_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis6GL>;

/// Q7 hexahedron (8^3 nodes).
using Q7_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis7GL>;

/// Q8 hexahedron (9^3 nodes).
using Q8_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis8GL>;

/// Q9 hexahedron (10^3 nodes).
using Q9_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis9GL>;

/**
 * @brief Maps a polynomial order to its hexahedron class, in member @c type.
 * @tparam ORDER Polynomial order, from 1 to 9; other values do not compile.
 */
template <int ORDER>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector;

/// @cond Doxygen_Suppress

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<1> {
  using type = Q1_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<2> {
  using type = Q2_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<3> {
  using type = Q3_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<4> {
  using type = Q4_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<5> {
  using type = Q5_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<6> {
  using type = Q6_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<7> {
  using type = Q7_Hexahedron_Lagrange_GaussLobatto;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<8> {
  using type = Q8_Hexahedron_Lagrange_GaussLobatto;
};
template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<9> {
  using type = Q9_Hexahedron_Lagrange_GaussLobatto;
};
/// @endcond
#if __GNUC__
#pragma GCC diagnostic pop
#endif
#undef PARENT_GRADIENT_METHOD

#endif  //_QkHEXAHEDRON_HPP_
