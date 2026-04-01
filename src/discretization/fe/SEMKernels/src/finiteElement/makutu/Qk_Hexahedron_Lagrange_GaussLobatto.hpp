/**
 * @file Qk_Hexahedron_Lagrange_GaussLobatto.hpp
 */

#ifndef _QkHEXAHEDRON_HPP_
#define _QkHEXAHEDRON_HPP_

#include <data_type.h>

#include "LagrangeBasis1.hpp"
#include "LagrangeBasis2.hpp"
#include "LagrangeBasis3GL.hpp"
#include "LagrangeBasis4GL.hpp"
#include "LagrangeBasis5GL.hpp"
#include "common/mathUtilites.hpp"

/**
 * This class is the basis class for the hexahedron finite element cells with
 * shape functions defined on Gauss-Lobatto quadrature points.
 * All the degree-specific versions (Q1, Q2, Q3, ...) are defined at the end of
 * this file.
 */
template <typename GL_BASIS>
class Qk_Hexahedron_Lagrange_GaussLobatto final
{
 public:
  constexpr static bool isShiva = false;

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
  constexpr static int numNodesPerFace =
      GL_BASIS::TensorProduct2D::numSupportPoints;

  /// The maximum number of support points per element.
  constexpr static int maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static int numQuadraturePoints = numNodes;

  struct TransformType
  {
    float data[8][3];
  };

  struct JacobianType
  {
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
  constexpr static int linearIndex3DVal(const int qa, int const qb,
                                        int const qc)
  {
    return qa + qb * num1dNodes + qc * numNodesPerFace;
  }

  /**
   * @brief Converts from the index of the point in the mesh and the linear 3D
   * index of the corresponding dof.
   * @param k The index of the mesh vertex, from 0 to 7
   * @return The linear index in 3D
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex3D(int const k)
  {
    return linearIndex3DVal((num1dNodes - 1) * (k % 2),
                            (num1dNodes - 1) * ((k % 4) / 2),
                            (num1dNodes - 1) * (k / 4));
  }

  /**
   * @brief The linear index associated to the given one-dimensional indices in
   * the two directions
   * @param qa The index in the first direction
   * @param qb The index in the second direction
   * @return The linear index in 2D
   */
  PROXY_HOST_DEVICE
  constexpr static int linearIndex2DVal(const int qa, const int qb)
  {
    return qa + qb * num1dNodes;
  }

  /**
   * @brief Converts from the index of the point in the mesh and the linear 2D
   * index of the corresponding dof.
   * @param k The index of the mesh vertex, from 0 to 3
   * @return The linear index in 2D
   */
  PROXY_HOST_DEVICE
  constexpr static int meshIndexToLinearIndex2D(int const k)
  {
    return linearIndex2DVal((num1dNodes - 1) * (k % 2),
                            (num1dNodes - 1) * (k / 2));
  }

  PROXY_HOST_DEVICE
  ~Qk_Hexahedron_Lagrange_GaussLobatto() = default;

  PROXY_HOST_DEVICE
  virtual int getNumQuadraturePoints()  // const override
  {
    return numQuadraturePoints;
  }

  PROXY_HOST_DEVICE
  virtual int getNumSupportPoints()  // const override
  {
    return numNodes;
  }

  PROXY_HOST_DEVICE
  virtual int getMaxSupportPoints() const { return maxSupportPoints; }

  /**
   * @brief Calculate shape functions values at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape
   * function value
   * @param[out] N The shape function values.
   */
  PROXY_HOST_DEVICE
  static void calcN(double const (&coords)[3], double (&N)[numNodes])
  {
    GL_BASIS::TensorProduct3D::value(coords, N);
  }

  /**
   * @brief Compute the interpolation coefficients of the q-th quadrature point
   * in a given direction
   * @param q the index of the quadrature point in 1D
   * @param k the index of the interval endpoint (0 or 1)
   * @return The interpolation coefficient
   */
  constexpr static real_t interpolationCoord(const int q, const int k)
  {
    const real_t alpha = (GL_BASIS::parentSupportCoord(q) + 1.0) / 2.0;
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
  constexpr static real_t basisGradientAt(const int q, const int p)
  {
    if (p <= halfNodes)
    {
      return GL_BASIS::gradientAt(q, p);
    }
    else
    {
      return -GL_BASIS::gradientAt(GL_BASIS::numSupportPoints - 1 - q,
                                   GL_BASIS::numSupportPoints - 1 - p);
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
  constexpr static real_t jacobianCoefficient1D(const int q, const int i,
                                                const int k, const int dir)
  {
    if (i == dir)
    {
      return k == 0 ? -1.0 / 2.0 : 1.0 / 2.0;
    }
    else
    {
      return interpolationCoord(q, k);
    }
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */

  PROXY_HOST_DEVICE
  static void calcN(int const q, real_t (&N)[numNodes])
  {
    for (int a = 0; a < numNodes; ++a)
    {
      N[a] = 0;
    }
    N[q] = 1.0;
  }

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the mesh support points.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */

  PROXY_HOST_DEVICE
  static real_t calcGradN(int const q, real_t const (&X)[numNodes][3],
                          real_t (&gradN)[numNodes][3]);
  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape
   * function value
   * @param[in] X Array containing the coordinates of the support points.
   * @param[out] gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */

  PROXY_HOST_DEVICE
  static real_t calcGradN(real_t const (&coords)[3],
                          real_t const (&X)[numNodes][3],
                          real_t (&gradN)[numNodes][3]);

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the mesh corners.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */

  PROXY_HOST_DEVICE
  static real_t calcGradNWithCorners(int const q, real_t const (&X)[8][3],
                                     real_t (&gradN)[numNodes][3]);
  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape
   * function value
   * @param[in] X Array containing the coordinates of the mesh corners.
   * @param[out] gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */

  PROXY_HOST_DEVICE
  static real_t calcGradNWithCorners(real_t const (&coords)[3],
                                     real_t const (&X)[8][3],
                                     real_t (&gradN)[numNodes][3]);

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  PROXY_HOST_DEVICE
  static real_t transformedQuadratureWeight(int const q,
                                            real_t const (&X)[numNodes][3]);

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space on a 2D domain
   * (face).
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation2d(int const qa, int const qb,
                                       real_t const (&X)[4][3],
                                       real_t (&J)[3][2]);

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianTransformation(int const qa, int const qb,
                                          int const qc, real_t const (&X)[8][3],
                                          real_t (&J)[3][3])
  {
    jacobianTransformation(qa, qb, qc, X, J);
    return invert3x3(J);
  }

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  PROXY_HOST_DEVICE
  static real_t invJacobianTransformation(int const q, real_t const (&X)[8][3],
                                          real_t (&J)[3][3])
  {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
    return invJacobianTransformation(qa, qb, qc, X, J);
  }

  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a quadrature point using the stored inverse of the Jacobian
   *   transformation matrix.
   * @param q The quadrature point index
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The symmetric gradient in Voigt notation.
   */
  PROXY_HOST_DEVICE
  static void symmetricGradient(int const q, real_t const (&invJ)[3][3],
                                real_t const (&var)[numNodes][3],
                                real_t (&grad)[6]);

  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param q The quadrature point index
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j}
   * var_{ai}\right ), \f]
   *
   */
  PROXY_HOST_DEVICE
  static void gradient(int const q, real_t const (&invJ)[3][3],
                       real_t const (&var)[numNodes][3], real_t (&grad)[3][3]);

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation(int const qa, int const qb, int const qc,
                                     real_t const (&X)[8][3],
                                     real_t (&J)[3][3]);

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space at a single
   * point.
   * @param coords The parent coordinates at which to evaluate the shape
   * function value
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformation(real_t const (&coords)[3],
                                     real_t const (&X)[numNodes][3],
                                     real_t (&J)[3][3]);

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space at a single
   * point. Assumes that the coordinate of high-order nodes are given by
   * trilinear interpolation of the mesh corners.
   * @param coords The parent coordinates at which to evaluate the shape
   * function value
   * @param X Array containing the coordinates of the mesh corners.
   * @param J Array to store the Jacobian transformation.
   */
  PROXY_HOST_DEVICE
  static void jacobianTransformationWithCorners(real_t const (&coords)[3],
                                                real_t const (&X)[8][3],
                                                real_t (&J)[3][3]);

  /**
   * @brief performs a trilinear interpolation to determine the real-world
   * coordinates of a vertex
   * @param[in] alpha Interpolation coefficient in [0,1] for the first
   * coordinate
   * @param[in] beta Interpolation coefficient in [0,1] for the second
   * coordinate
   * @param[in] gamma Interpolation coefficient in [0,1] for the third
   * coordinate
   * @param[in] X Real-world coordinates of the cell corners
   * @param[out] coords Real-world coordinates of the interpolated point
   */
  PROXY_HOST_DEVICE
  static void trilinearInterp(real_t const alpha, real_t const beta,
                              real_t const gamma, real_t const (&X)[8][3],
                              real_t (&coords)[3]);

  /**
   * @brief computes the real-world coordinates of the support nodes
   * @param[in] Xmesh Array containing the coordinates of the corners of the
   * mesh element
   * @param[out] X Array containing the coordinates of the support points.
   */
  PROXY_HOST_DEVICE
  static void computeLocalCoords(real_t const (&Xmesh)[8][3],
                                 real_t (&X)[numNodes][3]);

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   mass matrix M, i.e., the superposition matrix of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the mesh support points.
   * @return The diagonal mass term associated to q
   */
  template <typename FUNC>
  PROXY_HOST_DEVICE static void computeMassTerm(
      TransformType const &transformData, FUNC &&func);

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   damping matrix M, i.e., the superposition matrix of the shape functions
   *   integrated over a face.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @return The diagonal damping term associated to q
   */
  PROXY_HOST_DEVICE
  static real_t computeDampingTerm(int const q, real_t const (&X)[4][3]);

  /**
   * @brief computes the matrix B, defined as J^{-T}J^{-1}/det(J), where J is
   * the Jacobian matrix, at the given Gauss-Lobatto point.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian
   * @param B Array to store the matrix B, in Voigt notation
   */
  PROXY_HOST_DEVICE
  static void computeBMatrix(int const qa, int const qb, int const qc,
                             real_t const (&X)[8][3], real_t (&J)[3][3],
                             real_t (&B)[6]);

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexed by q to
   * the stiffness matrix R, i.e., the superposition matrix of first derivatives
   *   of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func1 Callback function accepting three parameters: qa,qb,qc invoked
   * when processing each quadrature point. This function will return each 1-D
   * quadrature used inside the kernel to compute the global degree of freedom
   * index, used when the model is defined on nodes.
   * @param func2 Callback function accepting three parameters: i, j and R_ij
   * invoked when processing each stiffness matrix contribution. The function
   *              will compute the index i (line of the matrix used to store the
   * result inside Y vector), index j (column of the matrix used to pick the
   * correct value from the output solution ) and R_ij which is the value of the
   * stiffness matrix itself.
   */
  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffnessTerm(
      TransformType const &transformData, FUNC1 &&func1, FUNC2 &&func2);

  /**
   * @brief Computes the "Grad(Phi)*B*Grad(Phi)" coefficient of the stiffness
   * term. The matrix B must be provided and Phi denotes a basis function.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param B Array of the B matrix, in Voigt notation
   * @param func1 Callback function accepting three parameters: qa,qb,qc invoked
   * when processing each quadrature point. This function will return each 1-D
   * quadrature used inside the kernel to compute the global degree of freedom
   * index, used when the model is defined on nodes.
   * @param func2 Callback function accepting three parameters: i, j and R_ij
   * invoked when processing each stiffness matrix contribution. The function
   *              will compute the index i (line of the matrix used to store the
   * result inside Y vector), index j (column of the matrix used to pick the
   * correct value from the output solution ) and R_ij which is the value of the
   * stiffness matrix itself.
   */
  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiBGradPhi(real_t const (&B)[6],
                                                       FUNC1 &&func1,
                                                       FUNC2 &&func2);

  /**
   * @brief Computes the "Grad(Phi)*Grad(Phi)" coefficient of the stiffness
   * term.
   * @tparam qa The 1D quadrature point index in xi0 direction (0,1)
   * @tparam qb The 1D quadrature point index in xi1 direction (0,1)
   * @tparam qc The 1D quadrature point index in xi2 direction (0,1)
   * @tparam FUNC1 First callback function type which takes four parameters:
   *               the three 1D Gauss-Lobatto point indices and the Jacobian
   * matrix
   * @tparam FUNC2 Second callback function type for processing computed
   * gradient products to get R_ij
   * @param[in] X Array of 8 nodal coordinates [node][dimension] defining the
   * corner of the hexaedra
   * @param[in,out] J Jacobian matrix [3][3] used for coordinate transformation
   * computations
   * @param[in] func1 First callback function invoked during gradient
   * computation
   * @param[in] func2 Second callback function invoked for gradient product
   * processing
   */
  template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeGradPhiGradPhi(JacobianType &J,
                                                      FUNC1 &&func1,
                                                      FUNC2 &&func2);

  /**
   * @brief Computes the non-zero contributions of the d.o.f. indexed by q to
   * the stiffness matrix R, i.e., the superposition matrix of first derivatives
   *        of the shape functions.
   * @tparam FUNC1 First callback function type invoked during Jacobian
   * computation
   * @tparam FUNC2 Second callback function type for processing stiffness
   * contributions
   * @param[in] X Array of 8 nodal coordinates [node][dimension] defining the
   * hexahedral element geometry
   * @param[in] func1 First callback function invoked for Jacobian-related
   * operations
   * @param[in] func2 Second callback function invoked to process computed
   * stiffness matrix contributions
   */
  template <typename FUNC1, typename FUNC2>
  PROXY_HOST_DEVICE static void computeStiffNessTermwithJac(
      TransformType const &transformData, FUNC1 &&func1, FUNC2 &&func2);

  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param q The quadrature point index
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  PROXY_HOST_DEVICE
  static void applyTransformationToParentGradients(
      int const q, real_t const (&invJ)[3][3], real_t (&gradN)[numNodes][3]);

  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space at a single point.
   * @param coords The parent coordinates at which to apply the transformation
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  PROXY_HOST_DEVICE
  static void applyTransformationToParentGradients(
      real_t const (&coords)[3], real_t const (&invJ)[3][3],
      real_t (&gradN)[numNodes][3]);

 private:
  /// The length of one dimension of the parent element.
  constexpr static real_t parentLength =
      GL_BASIS::parentSupportCoord(1) - GL_BASIS::parentSupportCoord(0);

  /// The volume of the element in the parent configuration.
  constexpr static real_t parentVolume =
      parentLength * parentLength * parentLength;
  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param coords The parent coordinates at which to evaluate the shape
   * function value
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template <typename FUNC, typename... PARAMS>
  PROXY_HOST_DEVICE static void supportLoop(real_t const (&coords)[3],
                                            FUNC &&func, PARAMS &&...params);
  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param q The quadrature node at which to evaluate the shape function value
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template <typename FUNC, typename... PARAMS>

  PROXY_HOST_DEVICE static void supportLoop(int const q, FUNC &&func,
                                            PARAMS &&...params);
};

/// @cond Doxygen_Suppress

template <typename GL_BASIS>
template <typename FUNC, typename... PARAMS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::supportLoop(
    real_t const (&coords)[3], FUNC &&func, PARAMS &&...params)
{
  for (int c = 0; c < num1dNodes; ++c)
  {
    for (int b = 0; b < num1dNodes; ++b)
    {
      for (int a = 0; a < num1dNodes; ++a)
      {
        real_t const dNdXi[3] = {
            GL_BASIS::gradient(a, coords[0]) * GL_BASIS::value(b, coords[1]) *
                GL_BASIS::value(c, coords[2]),
            GL_BASIS::value(a, coords[0]) * GL_BASIS::gradient(b, coords[1]) *
                GL_BASIS::value(c, coords[2]),
            GL_BASIS::value(a, coords[0]) * GL_BASIS::value(b, coords[1]) *
                GL_BASIS::gradient(c, coords[2])};

        int const nodeIndex = GL_BASIS::TensorProduct3D::linearIndex(a, b, c);

        func(dNdXi, nodeIndex, std::forward<PARAMS>(params)...);
      }
    }
  }
}

template <typename GL_BASIS>
template <typename FUNC, typename... PARAMS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::supportLoop(int const q,
                                                           FUNC &&func,
                                                           PARAMS &&...params)
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
  for (int c = 0; c < num1dNodes; ++c)
  {
    for (int b = 0; b < num1dNodes; ++b)
    {
      for (int a = 0; a < num1dNodes; ++a)
      {
        real_t const dNdXi[3] = {
            (b == qb && c == qc) ? basisGradientAt(a, qa) : 0,
            (a == qa && c == qc) ? basisGradientAt(b, qb) : 0,
            (a == qa && b == qb) ? basisGradientAt(c, qc) : 0};

        int const nodeIndex = GL_BASIS::TensorProduct3D::linearIndex(a, b, c);

        func(dNdXi, nodeIndex, std::forward<PARAMS>(params)...);
      }
    }
  }
}

//*************************************************************************************************

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradN(
    int const q, real_t const (&X)[numNodes][3], real_t (&gradN)[numNodes][3])
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
  real_t Xmesh[8][3] = {{0}};
  for (int k = 0; k < 8; k++)
  {
    const int nodeIndex = meshIndexToLinearIndex3D(k);
    for (int i = 0; i < 3; i++)
    {
      Xmesh[k][i] = X[nodeIndex][i];
    }
  }
  real_t J[3][3] = {{0}};

  jacobianTransformation(qa, qb, qc, Xmesh, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(q, J, gradN);

  return detJ;
}
//*************************************************************************************************
template <typename GL_BASIS>
PROXY_HOST_DEVICE

    real_t
    Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradN(
        real_t const (&coords)[3], real_t const (&X)[numNodes][3],
        real_t (&gradN)[numNodes][3])
{
  real_t J[3][3] = {{0}};

  jacobianTransformation(coords, X, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(coords, J, gradN);

  return detJ;
}
template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradNWithCorners(
    int const q, real_t const (&X)[8][3], real_t (&gradN)[numNodes][3])
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);

  real_t J[3][3] = {{0}};

  jacobianTransformation(qa, qb, qc, X, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(q, J, gradN);

  return detJ;
}
//*************************************************************************************************
template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::calcGradNWithCorners(
    real_t const (&coords)[3], real_t const (&X)[8][3],
    real_t (&gradN)[numNodes][3])
{
  real_t J[3][3] = {{0}};

  jacobianTransformationWithCorners(coords, X, J);

  real_t const detJ = invert3x3(J);

  applyTransformationToParentGradients(coords, J, gradN);

  return detJ;
}

//*************************************************************************************************
#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

// The following helper are used to compute a triple nested compile-time loop
// over the degrees of freedom in each direction of the hexahedral element
// when computing the stiffness and mass matrix contributions.

/*
 * Helper to perform a compile-time loop from 0 to N-1, calling a lambda with
 * std::integral_constant<int, I> as argument.
 */
template <int N, typename F, int... Is>
constexpr void for_constexpr_impl(F &&f, std::integer_sequence<int, Is...>)
{
  (f(std::integral_constant<int, Is>{}), ...);
}

/*
 * Perform a compile-time loop from 0 to N-1, calling a lambda with
 * std::integral_constant<int, I> as argument.
 */
template <int N, typename F>
constexpr void for_constexpr(F &&f)
{
  for_constexpr_impl<N>(std::forward<F>(f),
                        std::make_integer_sequence<int, N>{});
}

/*
 * Perform a triple nested compile-time loop from 0 to BoundI-1, 0 to BoundJ-1,
 * 0 to BoundK-1, calling a lambda with std::integral_constant<int, I>,
 * std::integral_constant<int, J>, std::integral_constant<int, K> as arguments.
 */

template <int BoundI, int BoundJ, int BoundK, typename Lambda>
constexpr void triple_loop(Lambda &&lambda)
{
  for_constexpr<BoundI>([&](auto I) {
    for_constexpr<BoundJ>([&](auto J) {
      for_constexpr<BoundK>([&](auto K) { lambda(I, J, K); });
    });
  });
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformation(
    int const qa, int const qb, int const qc, real_t const (&X)[8][3],
    real_t (&J)[3][3])
{
  for (int k = 0; k < 8; k++)
  {
    const int ka = k % 2;
    const int kb = (k % 4) / 2;
    const int kc = k / 4;
    for (int j = 0; j < 3; j++)
    {
      real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) *
                        jacobianCoefficient1D(qb, 1, kb, j) *
                        jacobianCoefficient1D(qc, 2, kc, j);
      for (int i = 0; i < 3; i++)
      {
        J[i][j] += jacCoeff * X[k][i];
      }
    }
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformation(
    real_t const (&coords)[3], real_t const (&X)[numNodes][3],
    real_t (&J)[3][3])
{
  supportLoop(
      coords,
      [](real_t const(&dNdXi)[3], int const nodeIndex,
         real_t const(&X)[numNodes][3], real_t(&J)[3][3]) {
        real_t const *Xnode = X[nodeIndex];
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            J[i][j] = J[i][j] + dNdXi[j] * Xnode[i];
          }
        }
      },
      X, J);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<
    GL_BASIS>::jacobianTransformationWithCorners(real_t const (&coords)[3],
                                                 real_t const (&X)[8][3],
                                                 real_t (&J)[3][3])
{
  supportLoop(
      coords,
      [](real_t const(&dNdXi)[3], int const nodeIndex, real_t const(&X)[8][3],
         real_t(&J)[3][3]) {
        int qa, qb, qc;
        GL_BASIS::TensorProduct3D::multiIndex(nodeIndex, qa, qb, qc);
        real_t Xnode[3];
        real_t alpha = (GL_BASIS::parentSupportCoord(qa) + 1.0) / 2.0;
        real_t beta = (GL_BASIS::parentSupportCoord(qb) + 1.0) / 2.0;
        real_t gamma = (GL_BASIS::parentSupportCoord(qc) + 1.0) / 2.0;
        trilinearInterp(alpha, beta, gamma, X, Xnode);
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            J[i][j] = J[i][j] + dNdXi[j] * Xnode[i];
          }
        }
      },
      X, J);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::trilinearInterp(
    real_t const alpha, real_t const beta, real_t const gamma,
    real_t const (&X)[8][3], real_t (&coords)[3])
{
  for (int i = 0; i < 3; i++)
  {
    coords[i] = X[0][i] * (1.0 - alpha) * (1.0 - beta) * (1.0 - gamma) +
                X[1][i] * alpha * (1.0 - beta) * (1.0 - gamma) +
                X[2][i] * (1.0 - alpha) * beta * (1.0 - gamma) +
                X[3][i] * alpha * beta * (1.0 - gamma) +
                X[4][i] * (1.0 - alpha) * (1.0 - beta) * gamma +
                X[5][i] * alpha * (1.0 - beta) * gamma +
                X[6][i] * (1.0 - alpha) * beta * gamma +
                X[7][i] * alpha * beta * gamma;
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeLocalCoords(
    real_t const (&Xmesh)[8][3], real_t (&X)[numNodes][3])
{
  int qa, qb, qc;
  for (int q = 0; q < numNodes; q++)
  {
    GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
    real_t alpha = (GL_BASIS::parentSupportCoord(qa) + 1.0) / 2.0;
    real_t beta = (GL_BASIS::parentSupportCoord(qb) + 1.0) / 2.0;
    real_t gamma = (GL_BASIS::parentSupportCoord(qc) + 1.0) / 2.0;
    trilinearInterp(alpha, beta, gamma, Xmesh, X[q]);
  }
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::jacobianTransformation2d(
    int const qa, int const qb, real_t const (&X)[4][3], real_t (&J)[3][2])
{
  for (int k = 0; k < 4; k++)
  {
    int ka = k % 2;
    int kb = k / 2;
    for (int j = 0; j < 2; j++)
    {
      real_t jacCoeff = jacobianCoefficient1D(qa, 0, ka, j) *
                        jacobianCoefficient1D(qb, 1, kb, j);
      for (int i = 0; i < 3; i++)
      {
        J[i][j] += jacCoeff * X[k][i];
      }
    }
  }
}

template <typename GL_BASIS>
template <typename FUNC>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeMassTerm(
    TransformType const &transformData, FUNC &&func)
{
  constexpr int N = num1dNodes;
  triple_loop<N, N, N>([&](auto const icqa, auto const icqb, auto const icqc) {
    constexpr int qa = decltype(icqa)::value;
    constexpr int qb = decltype(icqb)::value;
    constexpr int qc = decltype(icqc)::value;
    constexpr int q = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, qc);
    constexpr real_t w3D =
        GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);
    real_t J[3][3] = {{0}};
    float const(&X)[8][3] = transformData.data;
    jacobianTransformation(qa, qb, qc, X, J);
    real_t val = std::abs(determinant(J)) * w3D;
    func(q, val);
  });
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE real_t
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeDampingTerm(
    int const q, real_t const (&X)[4][3])
{
  int qa, qb;
  GL_BASIS::TensorProduct2D::multiIndex(q, qa, qb);
  const real_t w2D = GL_BASIS::weight(qa) * GL_BASIS::weight(qb);
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
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeBMatrix(
    int const qa, int const qb, int const qc, real_t const (&X)[8][3],
    real_t (&J)[3][3], real_t (&B)[6])
{
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
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeGradPhiBGradPhi(
    real_t const (&B)[6], FUNC1 &&func1, FUNC2 &&func2)
{
  constexpr real_t wa = GL_BASIS::weight(qa);
  constexpr real_t wb = GL_BASIS::weight(qb);
  constexpr real_t wc = GL_BASIS::weight(qc);
  const real_t w =
      GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);
  func1(qa, qb, qc);
  for (int i = 0; i < num1dNodes; i++)
  {
    const int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
    const int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
    const int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
    const real_t gia = basisGradientAt(i, qa);
    const real_t gib = basisGradientAt(i, qb);
    const real_t gic = basisGradientAt(i, qc);
    for (int j = 0; j < num1dNodes; j++)
    {
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
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeStiffnessTerm(
    TransformType const &transformData, FUNC1 &&func1, FUNC2 &&func2)
{
  triple_loop<num1dNodes, num1dNodes, num1dNodes>(
      [&](auto const icqa, auto const icqb, auto const icqc) {
        constexpr int qa = decltype(icqa)::value;
        constexpr int qb = decltype(icqb)::value;
        constexpr int qc = decltype(icqc)::value;
        real_t B[6] = {0};
        real_t J[3][3] = {{0}};
        float const(&X)[8][3] = transformData.data;
        computeBMatrix(qa, qb, qc, X, J, B);
        computeGradPhiBGradPhi<qa, qb, qc>(B, func1, func2);
      });
}

template <typename GL_BASIS>
template <typename FUNC1, typename FUNC2>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeStiffNessTermwithJac(
    TransformType const &transformData, FUNC1 &&func1, FUNC2 &&func2)
{
  triple_loop<num1dNodes, num1dNodes, num1dNodes>(
      [&](auto const icqa, auto const icqb, auto const icqc) {
        constexpr int qa = decltype(icqa)::value;
        constexpr int qb = decltype(icqb)::value;
        constexpr int qc = decltype(icqc)::value;
        JacobianType J = {{0}};
        jacobianTransformation(qa, qb, qc, transformData.data, J.data);
        computeGradPhiGradPhi<qa, qb, qc>(J, func1, func2);
      });
}

template <typename GL_BASIS>
template <int qa, int qb, int qc, typename FUNC1, typename FUNC2>
PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::computeGradPhiGradPhi(
    JacobianType &J, FUNC1 &&func1, FUNC2 &&func2)
{
  real_t const detJ = invert3x3(J.data);
  const real_t w =
      GL_BASIS::weight(qa) * GL_BASIS::weight(qb) * GL_BASIS::weight(qc);
  func1(qa, qb, qc, J.data);
  for (int i = 0; i < num1dNodes; i++)
  {
    const int ibc = GL_BASIS::TensorProduct3D::linearIndex(i, qb, qc);
    const int aic = GL_BASIS::TensorProduct3D::linearIndex(qa, i, qc);
    const int abi = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, i);
    const real_t gia = basisGradientAt(i, qa);
    const real_t gib = basisGradientAt(i, qb);
    const real_t gic = basisGradientAt(i, qc);
    for (int j = 0; j < num1dNodes; j++)
    {
      const int jbc = GL_BASIS::TensorProduct3D::linearIndex(j, qb, qc);
      const int ajc = GL_BASIS::TensorProduct3D::linearIndex(qa, j, qc);
      const int abj = GL_BASIS::TensorProduct3D::linearIndex(qa, qb, j);
      const real_t gja = basisGradientAt(j, qa);
      const real_t gjb = basisGradientAt(j, qb);
      const real_t gjc = basisGradientAt(j, qc);
      // diagonal terms
      const real_t w00 = w * gia * gja;
      func2(ibc, jbc, w00 * detJ, 0, 0);
      const real_t w11 = w * gib * gjb;
      func2(aic, ajc, w11 * detJ, 1, 1);
      const real_t w22 = w * gic * gjc;
      func2(abi, abj, w22 * detJ, 2, 2);
      // off-diagonal terms
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

//*************************************************************************************************
template <typename GL_BASIS>

PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::
    applyTransformationToParentGradients(int const q,
                                         real_t const (&invJ)[3][3],
                                         real_t (&gradN)[numNodes][3])
{
  supportLoop(
      q,
      [](real_t const(&dNdXi)[3], int const nodeIndex,
         real_t const(&invJ)[3][3], real_t(&gradN)[numNodes][3]) {
        // smaller register footprint by manually unrolling the for loops.
        gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] +
                              dNdXi[2] * invJ[2][0];
        gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] +
                              dNdXi[2] * invJ[2][1];
        gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] +
                              dNdXi[2] * invJ[2][2];
      },
      invJ, gradN);
}

//*************************************************************************************************
template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::
    applyTransformationToParentGradients(real_t const (&coords)[3],
                                         real_t const (&invJ)[3][3],
                                         real_t (&gradN)[numNodes][3])
{
  supportLoop(
      coords,
      [](real_t const(&dNdXi)[3], int const nodeIndex,
         real_t const(&invJ)[3][3], real_t(&gradN)[numNodes][3]) {
        gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] +
                              dNdXi[2] * invJ[2][0];
        gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] +
                              dNdXi[2] * invJ[2][1];
        gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] +
                              dNdXi[2] * invJ[2][2];
      },
      invJ, gradN);
}

template <typename GL_BASIS>

PROXY_HOST_DEVICE real_t
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::transformedQuadratureWeight(
    int const q, real_t const (&X)[numNodes][3])
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex(q, qa, qb, qc);
  real_t J[3][3] = {{0}};

  jacobianTransformation(qa, qb, qc, X, J);

  return determinant(J);
}

template <typename GL_BASIS>

PROXY_HOST_DEVICE void
Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::symmetricGradient(
    int const q, real_t const (&invJ)[3][3], real_t const (&var)[numNodes][3],
    real_t (&grad)[6])
{
  supportLoop(
      q,
      [](real_t const(&dNdXi)[3], int const nodeIndex,
         real_t const(&invJ)[3][3], real_t const(&var)[numNodes][3],
         real_t(&grad)[6]) {
        real_t gradN[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            gradN[i] = gradN[i] + dNdXi[j] * invJ[j][i];
          }
        }

        grad[0] = grad[0] + gradN[0] * var[nodeIndex][0];
        grad[1] = grad[1] + gradN[1] * var[nodeIndex][1];
        grad[2] = grad[2] + gradN[2] * var[nodeIndex][2];
        grad[3] = grad[3] + gradN[2] * var[nodeIndex][1] +
                  gradN[1] * var[nodeIndex][2];
        grad[4] = grad[4] + gradN[2] * var[nodeIndex][0] +
                  gradN[0] * var[nodeIndex][2];
        grad[5] = grad[5] + gradN[1] * var[nodeIndex][0] +
                  gradN[0] * var[nodeIndex][1];
      },
      invJ, var, grad);
}

template <typename GL_BASIS>
PROXY_HOST_DEVICE void Qk_Hexahedron_Lagrange_GaussLobatto<GL_BASIS>::gradient(
    int const q, real_t const (&invJ)[3][3], real_t const (&var)[numNodes][3],
    real_t (&grad)[3][3])
{
  supportLoop(
      q,
      [](real_t const(&dNdXi)[3], int const nodeIndex,
         real_t const(&invJ)[3][3], real_t const(&var)[numNodes][3],
         real_t(&grad)[3][3]) {
        for (int i = 0; i < 3; ++i)
        {
          real_t gradN = 0.0;
          ;
          for (int j = 0; j < 3; ++j)
          {
            gradN = gradN + dNdXi[j] * invJ[j][i];
          }
          for (int k = 0; k < 3; ++k)
          {
            grad[k][i] = grad[k][i] + gradN * var[nodeIndex][k];
          }
        }
      },
      invJ, var, grad);
}

//*************************************************************************************************
using Q1_Hexahedron_Lagrange_GaussLobatto =
    Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis1>;

using Q2_Hexahedron_Lagrange_GaussLobatto =
    Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis2>;

using Q3_Hexahedron_Lagrange_GaussLobatto =
    Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis3GL>;

using Q4_Hexahedron_Lagrange_GaussLobatto =
    Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis4GL>;

using Q5_Hexahedron_Lagrange_GaussLobatto =
    Qk_Hexahedron_Lagrange_GaussLobatto<LagrangeBasis5GL>;

template <int ORDER>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector;

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<1>
{
  using type = Q1_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<2>
{
  using type = Q2_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<3>
{
  using type = Q3_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<4>
{
  using type = Q4_Hexahedron_Lagrange_GaussLobatto;
};

template <>
struct Qk_Hexahedron_Lagrange_GaussLobatto_Selector<5>
{
  using type = Q5_Hexahedron_Lagrange_GaussLobatto;
};

#if __GNUC__
#pragma GCC diagnostic pop
#endif
#undef PARENT_GRADIENT_METHOD

#endif  //_QkHEXAHEDRON_HPP_
